#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:34:07 2020

@author: sw1906
"""
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tensorflow.python.framework import ops
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

# Set Seed
np.random.seed(0)

#Load in data
df = pd.read_csv('CNA_Full.csv', header = 0)
df = df.apply(pd.to_numeric, errors='ignore')
df = df.drop(df.columns[0], axis=1) # Remove Patient ID

# Use daatorg to speed up CNA_Full data analysis. 
# For CNA ER/PR/Intersection analyses, the CNA_ER/PR data was loaded in, and 
# the data was trained in the same manner, without the dataorg function. 

def dataorg(yval, specgroup = False):
    global df
    global train
    global test
    global Xtrain
    global Xtest
    global ytrain
    global ytest
    
    if yval == 'ER' or yval == 'PR':
        df = df.drop(df.columns[474:481], axis=1) #Remove Clin data
        
        if yval =='ER':
            if specgroup:
                df = df.loc[df['ER_Range'].isin(specgroup)]
            df = df.drop(df.columns[-2], axis=1) #Remove PR
            dfsplit = np.random.rand(len(df)) < 0.8
            train = df[dfsplit]
            test = df[~dfsplit]
            Xtrain = train.drop("ER_Range", axis=1)
            Xtest = test.drop("ER_Range", axis=1)
            ytrain = pd.get_dummies(train.ER_Range)
            ytest = pd.get_dummies(test.ER_Range)
        else: #PR has been chosen
            if specgroup:
                df = df.loc[df['PR_Range'].isin(specgroup)]
            df = df.drop(df.columns[-1], axis=1) #Remove ER
            dfsplit = np.random.rand(len(df)) < 0.8
            train = df[dfsplit]
            test = df[~dfsplit]
            Xtrain = train.drop("PR_Range", axis=1)
            Xtest = test.drop("PR_Range", axis=1)
            ytrain = pd.get_dummies(train.PR_Range)
            ytest = pd.get_dummies(test.PR_Range)
    elif yval == 'HRS':
        if specgroup:
            df = df.loc[df['HRS'].isin(specgroup)]
        df = df.drop(df.columns[475:481], axis=1) #Remove clin data
        df = df.drop(df.columns[-1], axis=1) #Remove PR
        df = df.drop(df.columns[-1], axis=1) #Remove ER
        dfsplit = np.random.rand(len(df)) < 0.8
        train = df[dfsplit]
        test = df[~dfsplit]
        Xtrain = train.drop("HRS", axis=1)
        Xtest = test.drop("HRS", axis=1)
        ytrain = pd.get_dummies(train.HRS)
        ytest = pd.get_dummies(test.HRS)
    else:
        sys.exit("Set a valid Y value. \n Use: 'ER', 'PR', or 'HRS' ")


# Create and train a tensorflow model of a neural network
def create_train_model(hidden_nodes, num_iters):
    
    # Reset the graph
    #tf.reset_default_graph()
    ops.reset_default_graph()

    # Placeholders for input and output data
    X = tf.placeholder(shape=(None, x_size), dtype=tf.float64, name='X')
    y = tf.placeholder(shape=(None, y_size), dtype=tf.float64, name='y')

    # Variables for two group of weights between the three layers of the network
    W1 = tf.Variable(np.random.rand(x_size, hidden_nodes), dtype=tf.float64)
    W2 = tf.Variable(np.random.rand(hidden_nodes, y_size), dtype=tf.float64)

    # Create the neural net graph
    A1 = tf.sigmoid(tf.matmul(X, W1))
    y_est = tf.sigmoid(tf.matmul(A1, W2))

    # Define a loss function
    deltas = tf.square(y_est - y)
    loss = tf.reduce_sum(deltas)

    # Define a train operation to minimize the loss
    optimizer = tf.train.GradientDescentOptimizer(0.005)
    train = optimizer.minimize(loss)

    # Initialize variables and run session
    init = tf.global_variables_initializer()
    sess = tf.Session()
    sess.run(init)

    # Go through num_iters iterations
    for i in range(num_iters):
        sess.run(train, feed_dict={X: Xtrain, y: ytrain})
        loss_plot[hidden_nodes].append(sess.run(loss, 
                feed_dict={X: Xtrain.values, y: ytrain.values}))
        weights1 = sess.run(W1)
        weights2 = sess.run(W2)
        
    print("loss (hidden nodes: %d, iterations: %d): %.2f" % (hidden_nodes, num_iters, loss_plot[hidden_nodes][-1]))
    sess.close()
    return weights1, weights2


# Choose data in use, train the networks 
    
# Define which y value is being predicted          
specificHRS = ['HR+/HER2-', 'Triple Negative', 'HR+/HER2+']
highlow = ['0-25','75-100']
dataorg('HRS')



# Define size for network description
x_size = len(Xtrain.columns)
y_size = len(ytrain.columns) 

# Plot the loss function over iterations
num_hidden_nodes = [5, 10, 20, 50]
loss_plot = {5: [], 10: [], 20: [], 50: []}
weights1 = {5: None, 10: None, 20: None, 50: None}
weights2 = {5: None, 10: None, 20: None,50: None}
num_iters = 2000

plt.figure(figsize=(12,8))
for hidden_nodes in num_hidden_nodes:
    weights1[hidden_nodes], weights2[hidden_nodes] = create_train_model(hidden_nodes, num_iters)
    plt.plot(range(num_iters), loss_plot[hidden_nodes], label="nn: x-%d-y" % hidden_nodes)
    
plt.xlabel('Iteration', fontsize=12)
plt.ylabel('Loss', fontsize=12)
plt.legend(fontsize=12)


# Evaluate models on the test set
X = tf.placeholder(shape=(None, x_size), dtype=tf.float64, name='X')
y = tf.placeholder(shape=(None, y_size), dtype=tf.float64, name='y')

for hidden_nodes in num_hidden_nodes:

    # Forward propagation
    W1 = tf.Variable(weights1[hidden_nodes])
    W2 = tf.Variable(weights2[hidden_nodes])
    A1 = tf.sigmoid(tf.matmul(X, W1))
    y_est = tf.sigmoid(tf.matmul(A1, W2))

    # Calculate the predicted outputs
    init = tf.global_variables_initializer()
    with tf.Session() as sess:
        sess.run(init)
        y_est_np = sess.run(y_est, feed_dict={X: Xtest, y: ytest})

    # Calculate the prediction accuracy
    correct = [estimate.argmax(axis=0) == target.argmax(axis=0) 
               for estimate, target in zip(y_est_np, ytest.values)]
    accuracy = 100 * sum(correct) / len(correct)
    print('Network architecture %d-%d-%d, accuracy: %.2f%%' % (x_size, hidden_nodes, y_size, accuracy))
    
    
    
    
    
    
    

