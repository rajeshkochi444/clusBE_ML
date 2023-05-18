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
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import wandb
from wandb.keras import WandbCallback

data_x = np.load('../DataXY/x_all_adsorb_firstrowTM.npy')
data_y = np.load('../DataXY/y_all_adsorb_firstrowTM.npy')


#print(data)
print(data_x.shape)
#print(len(data.shape))
print(data_y.shape)

data_x_reshape = data_x.reshape(data_x.shape[0], data_x.shape[1]*data_x.shape[2])
#print(data_reshape)
print(data_x_reshape.shape)

X = copy.deepcopy(data_x_reshape)
y = copy.deepcopy(data_y)


# Select the hyperparameters you want to tune. This is specified like the following:
sweep_config = {
    #'method': 'random',
    'method': 'grid',
    #'name': 'NN sweep',
    #'metric': {
        #'goal': 'minimize',
        #'name': 'validation_loss'
	#	},
    'parameters': {
        'batch_size': {'values': [256]},
        #'epochs': {'values': [5, 10, 15]},
        'lr': {'values': [ 0.0001]},
        'layer1': {'values': [ 2048, 1024]},
        'layer2': {'values': [ 1024, 512, 256]},
        'layer3': {'values': [ 1024, 512, 256, 128]},
        'drop_rate': {'values': [  0.3, 0.4]},
        #'drop_rate': {'values': [ 0.2, 0.3, 0.4, 0.5, 0.6]},
     }
}

#initialize your sweep
#sweep_id = wandb.sweep(sweep=sweep_config,  project = 'clusnet')
sweep_id = wandb.sweep(sweep=sweep_config)

def train():
    wandb.init( project='clusnet')


    config = wandb.config # Config is a variable that holds and saves hyperparameters and inputs
    #config.learning_rate = 0.001
    config.epochs = 200
    #config.batch_size = 128
    config.random_seed = 88

    #config.layer1 = 1024
    #config.layer2 = 512
    
    X_train, X_val, y_train, y_val = train_test_split(X, y,  test_size=0.20, random_state=config.random_seed)
    
    input_shape = X_train.shape[1]

    model = Sequential()
    model.add(Dense(config.layer1, activation = 'relu', input_shape=(input_shape,), kernel_regularizer=keras.regularizers.l2(0.01)))
    model.add(Dropout(config.drop_rate))
    model.add(Dense(config.layer2, activation = 'relu', kernel_regularizer=keras.regularizers.l2(0.01)))
    model.add(Dropout(config.drop_rate))
    model.add(Dense(config.layer3, activation = 'relu', kernel_regularizer=keras.regularizers.l2(0.01)))
    model.add(Dropout(config.drop_rate))
    model.add(Dense(1, activation = 'linear', kernel_regularizer=keras.regularizers.l2(0.01)))

    model.summary()

    opt = keras.optimizers.Adam(lr=config.lr)

    #model.compile(optimizer=opt, loss="mse", metrics=['mse', 'mae'])
    model.compile(optimizer=opt, loss="mse" )


    model.fit(X_train, y_train,
                    epochs = config.epochs,
                    batch_size = config.batch_size,
                    validation_data = (X_val, y_val),
                    shuffle = True,
                    callbacks=[WandbCallback(validation_data=(X_val, y_val)), tf.keras.callbacks.EarlyStopping(patience=10, restore_best_weights=True)])

    #wandb.log({'val_loss': validation_loss 
     #        }
      #      )

wandb.agent(sweep_id, function=train)
#wandb.agent(sweep_id, function=train, count=5)

