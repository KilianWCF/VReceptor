#!/usr/bin/env python
# coding: utf-8

# LSTM model with wrapper for peptide-activity regression tasks.
# This code defines a PyTorch LSTM model and a wrapper class for training and prediction. The model
# is designed to handle sequences of features, such as those found in peptide-activity datasets.
# It includes methods for fitting the model to training data and making predictions on new data.
# written by Kilian Conde-Frieboes, version below is from 2025-07-01


import torch
import torch.nn as nn
from torch.utils.data import random_split, TensorDataset, DataLoader
import time


class LSTM(nn.Module):
    '''LSTM model for sequence regression tasks.
    This model uses an LSTM layer followed by a fully connected layer to output predictions.
    Args:
        input_size (int): Number of features in the input sequence.
        output_size (int): Number of output features.
        hiddensize (int, optional): Size of the hidden state in the LSTM. Defaults to 50.
        layers (int, optional): Number of LSTM layers. Defaults to 2.
        bidirectional (bool, optional): Whether to use a bidirectional LSTM. Defaults to True.
        FFNN_dropout (float, optional): Dropout rate for the fully connected layer. Defaults to 0.0.
        LSTM_dropout (float, optional): Dropout rate for the LSTM layers. Defaults to 0.0.
    '''
    def __init__(self, input_size, output_size, hiddensize = 50, layers = 2, bidirectional = True, FFNN_dropout = 0.0, LSTM_dropout = 0.0):
        super(LSTM, self).__init__()
        # Initialize the LSTM model with the specified parameters.
        
        self.hiddensize = hiddensize
        self.layers = layers
        self.bidirectional = bidirectional
        
        self.lstm = nn.LSTM (input_size = input_size, hidden_size = hiddensize, num_layers = layers, bidirectional = bidirectional, dropout = LSTM_dropout)
        self.ffnn_drop = nn.Dropout(FFNN_dropout)
        #self.ffnn = nn.Linear(in_features = hiddensize * (bidirectional + 1) * layers, out_features = 100)
        #self.l_out = nn.Linear(in_features=hiddensize*2*layers,out_features=3)
        self.l_out = nn.Linear(in_features = hiddensize * (bidirectional + 1) * layers, out_features = output_size)
        #nn.init.xavier_uniform_(self.ffnn.weight)
        
        nn.init.xavier_uniform_(self.l_out.weight)
       
    def forward(self, x):
        # Forward pass through the LSTM model.
        x = x.transpose(0,1)
        #input shape is (seq length, batch size, features)
        _, (h_n, c_n) = self.lstm(x)
        out = h_n.transpose(0, 1)
        out = out.reshape(-1, self.hiddensize * (self.bidirectional + 1) * self.layers)
        #out = self.ffnn(out)
        
        #out = torch.tanh(out)
        #out=self.ffnn_drop(out)
        out = self.l_out(out)
        return out
    
    def reset(self):
        self.lstm.reset_parameters()        
        self.l_out.reset_parameters()
        
        
        

    
    
class NN():
    '''Wrapper class for training and predicting with an LSTM model.
    This class handles the training process, including data loading, loss calculation, and optimization.
    It also provides methods for making predictions on new data.
    Args:
        model (LSTM): The LSTM model to be trained and used for predictions.
        GPU (bool, optional): Whether to use a GPU for training. Defaults to True.
    '''
    def __init__(self, model, GPU = True):
        self.GPU = GPU
        self.model=model
        
        if GPU:
            self.model.cuda()

    def fit(self, X, y, lr = 0.001, weight_decay = 0.001, max_iter = 1000, Verbose = True, reporting = 100, batch_size = 100):
        '''Fit the LSTM model to the training data.
        Args:
            X (array-like): Input features for training, shape (n_samples, n_features).
            y (array-like): Target values for training, shape (n_samples, n_outputs).
            lr (float, optional): Learning rate for the optimizer. Defaults to 0.001.
            weight_decay (float, optional): Weight decay for the optimizer. Defaults to 0.001.
            max_iter (int, optional): Maximum number of training iterations. Defaults to 1000.
            Verbose (bool, optional): Whether to print training progress. Defaults to True.
            reporting (int, optional): Frequency of reporting training progress. Defaults to 100.
            batch_size (int, optional): Batch size for training. Defaults to 100.
        '''
        self.model.train()
        X = torch.tensor(X, dtype=torch.float)
        y = torch.tensor(y, dtype=torch.float)
        
        train=TensorDataset(X,y)
        train_loader = DataLoader(train, batch_size = batch_size, shuffle = False,  drop_last = False)
        criterion = nn.MSELoss()
        optimizer = torch.optim.Adam(self.model.parameters(), lr = lr, weight_decay = weight_decay)

        start_time=time.time()
        
        for epoch in range(max_iter):
            for X,y in train_loader:
                if self.GPU:
                    X = X.cuda()
                    y = y.cuda()
                y_pre = self.model(X)
                train_loss = criterion(y_pre, y)
                optimizer.zero_grad()
                train_loss.backward()
                optimizer.step()
            if ((epoch+1)%reporting == 0) and (Verbose == True):
                print (f'Epoch: {epoch+1}\tTime since start: {time.time()-start_time:10.2f}s\tTrain loss: {train_loss:.3f}')
        
        
    def predict(self, X):
        '''Make predictions using the trained LSTM model.
        Args:
            X (array-like): Input features for prediction, shape (n_samples, n_features).
        Returns:
            numpy.ndarray: Predicted values, shape (n_samples, n_outputs).
        '''
        self.model.eval()
        X = torch.tensor(X, dtype=torch.float)
        if self.GPU:
            X=X.cuda()
        pre = self.model(X)
        self.model.train()
        return pre.cpu().detach().numpy()