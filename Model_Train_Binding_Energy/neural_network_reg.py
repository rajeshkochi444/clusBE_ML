import ase
import os
import glob
import numpy as np
from ase.io import read, write, Trajectory
import seaborn as sns
import matplotlib.pyplot as plt
import random
import itertools
import copy
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
import tensorflow as tf
from datetime import datetime
from sklearn.metrics import mean_squared_error, mean_absolute_error, explained_variance_score
import matplotlib.pyplot as plt
import wandb
from wandb.keras import WandbCallback

#config = dict(
        #architecture = "NN",
        #dataset_id = "filtered2eV",
        #)

wandb.init(
        project='clusnet_random3k',
        notes='1stRowTMs',
        tags=['1atomdes','random', '3k' 'NN'],
        )
dt_string = datetime.now().strftime("%d-%m-%Y-%H-%M-%S")
#wandb.run.name = dt_string
#wandb.run.save()
# Default values for hyper-parameters
config = wandb.config # Config is a variable that holds and saves hyperparameters and inputs
config.learning_rate = 0.001
config.epochs = 200
config.batch_size = 128
config.architecture_name = 'NN'
config.architecture = 1024,512,256,64
config.dropout  = 0.4, 0.4, 0.4, 0.4
config.early_stopping = 'Yes'
config.dataset_id = "acsf_1atom"
config.random_seed = 88
config.datasplit = 80,20
config.regularizer = 'L2'

wandb.run.name = config.architecture_name + '_' + config.regularizer + '_' +  config.dataset_id + '_seed' +  str(config.random_seed) + '_' + dt_string
wandb.run.save()

data_x = np.load('../../DataXY_Random/x_all_adsorb_firstrowTM.npy')
data_y = np.load('../../DataXY_Random/y_all_adsorb_firstrowTM.npy')

#print(data)
print(data_x.shape)
#print(len(data.shape))
print(data_y.shape)

data_x_reshape = data_x.reshape(data_x.shape[0], data_x.shape[1]*data_x.shape[2])
#print(data_reshape)
print(data_x_reshape.shape)

X = copy.deepcopy(data_x_reshape)
y = copy.deepcopy(data_y)

from sklearn.model_selection import train_test_split
X_train, X_val, y_train, y_val = train_test_split(X, y,  test_size=0.20, random_state=config.random_seed)

print(X_train.shape)
print(X_val.shape)
print(y_train.shape)
print(y_val.shape)


input_shape = X_train.shape[1]

model = Sequential()
model.add(Dense(1024, activation = 'relu', input_shape=(input_shape,), kernel_regularizer=keras.regularizers.l2(0.01)))
model.add(Dropout(0.4))
#model.add(Dense(512, activation = 'relu', input_dim=X_train.shape[1] ))
model.add(Dense(512, activation = 'relu', kernel_regularizer=keras.regularizers.l2(0.01) ))
model.add(Dropout(0.4))
model.add(Dense(256, activation = 'relu', kernel_regularizer=keras.regularizers.l2(0.01) ))
model.add(Dropout(0.4))
model.add(Dense(128, activation = 'relu', kernel_regularizer=keras.regularizers.l2(0.01) ))
model.add(Dropout(0.4))
model.add(Dense(64, activation = 'relu', kernel_regularizer=keras.regularizers.l2(0.01)))
model.add(Dropout(0.4))
model.add(Dense(1, activation = 'linear', kernel_regularizer=keras.regularizers.l2(0.01)))

model.summary()

opt = keras.optimizers.Adam(learning_rate=config.learning_rate)

model.compile(optimizer=opt, loss="mse", metrics=['mse', 'mae'])


#es = keras.callbacks.EarlyStopping(
#        monitor="val_loss", # metrics to monitor
#        patience=10, # how many epochs before stop
#        verbose=1,
#        restore_best_weights=True)

history = model.fit(X_train, y_train,
                    epochs = config.epochs,
                    batch_size = config.batch_size,
                    validation_data = (X_val, y_val),
                    shuffle = True,
                    callbacks=[WandbCallback(validation_data=(X_val, y_val)), tf.keras.callbacks.EarlyStopping(patience=10, restore_best_weights=True)])

predictions = model.predict(X_val)
mean_sq_error = mean_squared_error(y_val, predictions)
mean_abs_error = mean_absolute_error(y_val, predictions)
mean_rms_error = np.sqrt(mean_sq_error) 
explained_variance = explained_variance_score(y_val, predictions)
wandb.log({'val_mean_squared_error': mean_sq_error, 
           'val_mean_abs_error': mean_abs_error,
           'val_rms_error:': mean_rms_error,
           'explained_variance': explained_variance,
          }
         )
print('before plot')
print(y_val.shape)
print(predictions.shape)

y_val_flatten = y_val.flatten()
y_predict_flatten = predictions.flatten()
print(y_val_flatten.shape)
print(y_predict_flatten.shape)

data = [[x, y] for (x, y) in zip(y_val_flatten, y_predict_flatten)]
table = wandb.Table(data=data, columns = ["y_true", "y_predict"])
wandb.log({"scatter_y_val_predict" : wandb.plot.scatter(table,  "y_true", "y_predict")})

plt.figure(figsize=(12,6))
plt.scatter(y_val_flatten, y_predict_flatten)
plt.plot(y_val_flatten, y_val_flatten, 'r')
plt.title('y_val vs y_predict')
plt.xlabel('y_val')
plt.ylabel('y_predict')
plt.savefig('y_scatter.png')
wandb.log({"y_scatter_plot": wandb.Image("y_scatter.png")})

