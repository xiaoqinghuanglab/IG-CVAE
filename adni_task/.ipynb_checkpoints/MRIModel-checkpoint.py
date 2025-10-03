import os
import nibabel as nib
import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from scipy.ndimage import zoom
from numba import jit, cuda

class MRIDataset(Dataset):
    def __init__(self, annotations, root_dir, target_shape=(128, 128, 128),transform=None):
        self.annotations = annotations
        self.root_dir = root_dir
        self.transform = transform
        self.label_map = {'CN': 0, 'AD': 1}
        self.target_shape = target_shape

    def __len__(self):
        return len(self.annotations)

    def __getitem__(self, idx):
        image_id = self.annotations.iloc[idx]['Image Data ID']
        label = self.annotations.iloc[idx]['Group']

        img_path = self._find_nii_file(image_id)
        if not img_path:
            raise FileNotFoundError(f'NIfTI file not found for image_id: {image_id}')

        # image = nib.load(img_path).get_fdata()

        try:
            image = nib.load(img_path).get_fdata()
        except Exception as e:
            print(f'Error loading {img_path}: {e}')
            return self.__getitem__((idx + 1) % len(self))

        # Resize the image to the target shape
        try:
            image = self._resize_image(image, self.target_shape)
        except ValueError as e:
            print(f'Error resizing image {img_path}: {e}')
            # self.corrupted_files.append(img_path)
            return self.__getitem__((idx + 1) % len(self))

        if self.transform:
            image = self.transform(image)

        image = torch.tensor(image, dtype=torch.float32).unsqueeze(0)  # Add channel dimension
        label = self.label_map[label]
        label = torch.tensor(label, dtype=torch.long)

        return image, label

    def _find_nii_file(self, image_id):
        for root, _, files in os.walk(self.root_dir):
            if image_id in root:
                for file in files:
                    if file.endswith('.nii'):
                        return os.path.join(root, file)
        return None

    def _resize_image(self, image, target_shape):
        if image.ndim != len(target_shape):
            raise ValueError(f"Image dimensions {image.ndim} do not match target shape dimensions {len(target_shape)}")
        factors = [t / s for t, s in zip(target_shape, image.shape)]
        resized_image = zoom(image, factors, order=1)  # Use order=1 (bilinear interpolation)
        return resized_image

# load dataset
df = pd.read_csv("CSVData/Merged_MRI_Data.csv")
filtered_df = df[(df['Gene_Flag'] == 1) & (df['MRI_Flag'] == 1)]

# Remove duplicates based on the image ID
# filtered_df = filtered_df.drop_duplicates(subset='PTID')
filtered_df = filtered_df.dropna(subset=['DX'])
print(filtered_df['Image Data ID'].head(), filtered_df['PTID'].head(), filtered_df['DX'], filtered_df['Group'].head())
count_by_phase = filtered_df[(filtered_df['Gene_Flag'] == 1) & (filtered_df['MRI_Flag'] == 1)].groupby(
        'ORIGPROT').size().reset_index(name='Count')
print("\nCount by phase where both Gene_Flag and MRI_Flag are 1:")
print(count_by_phase)
# exit()
# Load the dataset
dataset = MRIDataset(filtered_df, root_dir='ADNI')
# print(dataset)

# Split dataset into train and test
train_size = int(0.8 * len(dataset))
test_size = len(dataset) - train_size
train_dataset, test_dataset = torch.utils.data.random_split(dataset, [train_size, test_size])

# Create data loaders
train_loader = DataLoader(train_dataset, batch_size=4, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=4, shuffle=False)


class CNN3DModel(nn.Module):
    def __init__(self):
        super(CNN3DModel, self).__init__()
        self.conv1 = nn.Conv3d(1, 32, kernel_size=3, stride=1, padding=1)
        self.conv2 = nn.Conv3d(32, 64, kernel_size=3, stride=1, padding=1)
        self.conv3 = nn.Conv3d(64, 128, kernel_size=3, stride=1, padding=1)
        self.pool = nn.MaxPool3d(kernel_size=2, stride=2, padding=0)
        self.dropout = nn.Dropout(p=0.5)

        # Calculate the size of the flattened features
        # Assuming input shape is (1, 128, 128, 128)
        input_shape = (1, 128, 128, 128)
        dummy_input = torch.zeros(1, *input_shape)
        dummy_output = self._forward_features(dummy_input)
        self.flattened_size = dummy_output.view(dummy_output.size(0), -1).size(1)

        self.fc1 = nn.Linear(self.flattened_size, 256)
        self.fc2 = nn.Linear(256, 2)

    def _forward_features(self, x):
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        x = self.pool(F.relu(self.conv3(x)))
        return x

    def forward(self, x):
        x = self._forward_features(x)
        x = x.view(x.size(0), -1)  # Flatten the features
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x)
        return x
def runMode():
    # Initialize the model, loss function, and optimizer
    model = CNN3DModel()
    model = nn.DataParallel(model)

    # Move the model to GPU
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)

    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    # Training loop
    num_epochs = 100
    print(train_loader)
    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        for inputs, labels in train_loader:
            # print(inputs, labels) jan cn - nov ad
            inputs, labels = inputs.to(device), labels.to(device)
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {running_loss/len(train_loader):.4f}')

    # Evaluation
    model.eval()
    correct = 0
    total = 0
    with torch.no_grad():
        for inputs, labels in test_loader:
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)
            _, predicted = torch.max(outputs.data, 1)
            total += labels.size(0)
            correct += (predicted == labels).sum().item()

    print(f'Accuracy: {100 * correct / total:.2f}%')
runMode()