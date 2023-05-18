import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import copy
from sklearn.metrics import mean_squared_error, mean_absolute_error, explained_variance_score
import matplotlib.pyplot as plt
from sklearn.preprocessing import  MinMaxScaler, StandardScaler
from sklearn.decomposition import PCA


data_x = np.load('acsf_ndes1_all_adsorb_1eV_Cu.npy')
data_y = np.load('y_all_adsorb.npy')


print(data_x.shape, data_y.shape)

data_x_reshape = data_x.reshape(data_x.shape[0], data_x.shape[1]*data_x.shape[2])
print(data_x_reshape.shape)


X = copy.deepcopy(data_x_reshape)
y = copy.deepcopy(data_y)

from sklearn.model_selection import train_test_split
X_train, X_val, y_train, y_val = train_test_split(X, y,  test_size=0.20, random_state=88)

print(X_train.shape)
print(X_val.shape)
print(y_train.shape)
print(y_val.shape)

#Normalize the data
scaler = StandardScaler()
X_train  = scaler.fit_transform(X_train)
X_val  = scaler.transform(X_val)

#PCA
pca = PCA(n_components=0.999)
pca.fit(X_train)


X_train_reduced = pca.transform(X_train)
X_val_reduced = pca.transform(X_val)

np.save('X_train_reduced.npy', X_train_reduced)
np.save('X_val_reduced.npy', X_val_reduced)
np.save('y_train.npy', y_train)
np.save('y_val.npy', y_val)


plt.grid()
plt.plot(np.cumsum(pca.explained_variance_ratio_ * 100))
plt.xlabel('Number of components')
plt.ylabel('Explained variance')
plt.savefig('cumsum_var_pca.png')

print(pca.explained_variance_ratio_ * 100)
print(np.cumsum(pca.explained_variance_ratio_ * 100))
cumsum_var = np.cumsum(pca.explained_variance_ratio_ * 100).tolist()
#print(cumsum_var)

for i, value in enumerate(cumsum_var):
    print(i, value)
