import warnings
import math
import numpy as np
import pylab
import matplotlib.pyplot as plt
from sklearn import linear_model
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error
from pandas import DataFrame
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.linear_model import LogisticRegression
import pandas as pd
from sklearn.metrics import r2_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
# Import function to automatically create polynomial features!
from sklearn.preprocessing import PolynomialFeatures
# Import Linear Regression and a regularized regression function
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LassoCV
from sklearn.pipeline import make_pipeline
from sklearn.tree import DecisionTreeRegressor

global finals, two, three, four, five, six,seven, eight, error, calc
finals, two, three, four, five, six, seven, eight, error, calc = [], [], [], [], [], [], [], [], [], []
with warnings.catch_warnings():
   warnings.simplefilter("ignore")
   from google.cloud import storage
   def main():
      '''tasks = ['bam2ta', 'bowtie2', 'filter', 'idr_pr',
      'macs2', 'macs2_pr1', 'macs2_pr2', 'macs2_signal_track', 'overlap_pr',
      'qc_report', 'read_genome_tsv', 'reproducibility_idr', 'reproducibility_overlap', 'spr',
      'trim_adapter', 'xcor']'''

      tasks = ['ataqc', 'bam2ta', 'bowtie2', 'filter', 'idr_pr',
      'macs2', 'macs2_pr1', 'macs2_pr2', 'macs2_signal_track', 'overlap_pr',
      'qc_report', 'read_genome_tsv', 'reproducibility_idr', 'reproducibility_overlap', 'spr',
      'trim_adapter', 'xcor']

      workflow_ids = [
         '92ecee2e-e74f-46ba-8c73-ed062687e434',
         '18ac4af8-6a5c-4970-a4ec-1d5c7f9e9d53',
         '8cf23c7a-d2cd-4d88-9d2a-42314709a1af',
         '61a4857b-25e5-4b88-9118-b3ac2e4f5aa5',
         '420405aa-9128-403c-aaf8-0c8fc62229c0',
         'eefc0137-1978-4a70-8b98-71eb90875f7f'
         ]

      sizes = [
         2,
         5,
         10,
         15,
         20,
         30
      ]
      client = storage.Client()
      bucket = client.get_bucket('gbsc-gcp-project-annohive-dev-user-dtluan-cromwell-test')
      result_dir = 'ENCSR356KRQ_subsampled_short/atac/'
      for task in tasks:
         # create data frame
         mems = []
         for i, _ in enumerate(sizes):
            mems.append(get_memory_usage(result_dir, workflow_ids[i], task, bucket))
         data = {
            'Size': sizes,
            'Memory': mems
         }
         plot_memory_usage(data, task)

   def plot_memory_usage(data, task):
      """ plot task memory usage of downsamples and predict memory usage of real sample
         Args:
               data: a table of sizes and memory usage arraya
               task: task name
         Returns: none
      """
      df = DataFrame(data, columns=['Size', 'Memory'])
      X = df[['Size']]
      Y = df['Memory']
      target = 100
      x_simulate = np.linspace(1, target, num=100)
      x_simulate = x_simulate.reshape(100, 1)
      ax = plt.subplot(1,1,1)
      plt.title(task)
      plt.xlabel('Down Sample Size (%)')
      plt.ylabel('Memory Usage (KB)')
      # add test data
      #plt.scatter(X, Y, color='red')
      # add linear regression
      score = linear_regression(plt, X, Y, x_simulate, target, ax)

      # add logarithmic and polynomial regression if linear regresion score < 0.9

      log_regression(plt, X, Y, x_simulate, target)
      poly_regression(plt, X, Y, x_simulate, target)
      # save figure
      # add test data
      plt.scatter(X, Y, color='red')

      #plot final point
      client = storage.Client()
      try:
          final_memory = get_memory_usage('ENCSR356KRQ_subsampled_short/atac/', 'a13a955a-1f6f-4376-a438-dc5dd9234efe', task, client.get_bucket('gbsc-gcp-project-annohive-dev-user-dtluan-cromwell-test'))
          finals.append(final_memory)
          plt.scatter(target, final_memory, color='pink')
      except:
          pass

      plt.savefig('plot/' + task + '.png')
      plt.clf()
   def linear_regression(plt, X, Y, x_simulate, target, ax):
      # create linear regression model
      regr = linear_model.LinearRegression()
      regr.fit(X, Y)
      yhat = regr.predict(X)
      y_simulate_linear = regr.predict(x_simulate)
      #plt.plot(X, yhat, color='blue')
      plt.plot(x_simulate, y_simulate_linear, color='blue', label = 'Linear')
      # predict target memory usage. input should be a 2D array
      target_memory = regr.predict([[ target ]])
      plt.scatter(target, target_memory[0], color='blue')
      plt.legend()
      #prediction = 'real sample: ' + str(int(target_memory[0]))  + ' MB'
      # get r-squared
      #score = regr.score(X, Y)
      #r_squared =  '$R^2$ = ' + str('{:0.5f}'.format(score))
      # get mean squared error
      #mse = mean_squared_error(yhat, Y)
      #rmse = 'RMSE = ' + str(int(math.sqrt(mse))) + ' MB'
      # get equation
      #coefficient = regr.coef_
      #intercept = regr.intercept_
      #equation = 'Y = ' + str(int(intercept)) + ' + ' + str(coefficient[0]) + 'X'
      #plt.text(0.5, 0.9, r_squared + '\n' + rmse + '\n' + prediction, ha='center', va='center', transform = ax.transAxes)
      #return score
   def poly_regression(plt, X, Y, x_simulate, target):
      y_formatted = np.array(Y)
      X_prediction = np.array([[100]])

      #lr = RandomForestClassifier()
      #lr.fit(X, y_formatted)
      #X_prediction = np.array([[100]])
      #predictions = lr.predict(X_prediction)
      #y_simulate_poly = lr.predict(x_simulate)
      #plt.plot(x_simulate, y_simulate_poly, color='orange')
      #plt.scatter(100, predictions, color = 'orange')
      #print("RandomForestClassifier error is ")
      #print(predictions - final_memory/final_memory*100)'''

     #lr = LogisticRegression()
      #lr.fit(X, y_formatted)
      #X_prediction = np.array([[100]])
      #predictions = lr.predict(X_prediction)
      #y_simulate_poly = lr.predict(x_simulate)
      #plt.plot(x_simulate, y_simulate_poly, color='orange')
      #plt.scatter(100, predictions, color = 'orange')
      #print("LogisticRegression error is ")
      #print(predictions - final_memory/final_memory*100)'''

      #lr = SVC()
      #lr.fit(X, y_formatted)
      #X_prediction = np.array([[100]])
      #predictions = lr.predict(X_prediction)
      #y_simulate_poly = lr.predict(x_simulate)
      #plt.plot(x_simulate, y_simulate_poly, color='orange')
      #plt.scatter(100, predictions, color = 'orange')'''


      # Alpha (regularization strength) of LASSO regression
      lasso_eps = 0.0001
      lasso_nalpha=20
      lasso_iter=5000
      # Min and max degree of polynomials features to consider
      degree_min = 2
      degree_max = 8
      y_formatted = np.array(Y)
      X_prediction = np.array([[100]])


      '''for degree in range(degree_min,degree_max+1):
          model = make_pipeline(PolynomialFeatures(degree, interaction_only=False), LassoCV(eps=lasso_eps,n_alphas=lasso_nalpha,max_iter=lasso_iter,
          normalize=True,cv=5))
          model.fit(X,y_formatted)
          X_prediction = np.array([[100]])
          test_pred2 = np.array(model.predict(X_prediction))
          test_pred = model.predict(x_simulate)
          if(degree == 2):
              two.append(test_pred)
              plt.scatter(100, test_pred2, color='magenta')
              plt.plot(x_simulate, test_pred, color='magenta' ,label = 'Degree 2')
              plt.legend()
          elif(degree == 3):
              three.append(test_pred)
              plt.scatter(100, test_pred2, color='green')
              plt.plot(x_simulate, test_pred, color='green', label = 'Degree 3')
              plt.legend()
          elif(degree == 4):
              four.append(test_pred)
              plt.scatter(100, test_pred2, color='blue')
              plt.plot(x_simulate, test_pred, color='blue', label = 'Degree 4')
              plt.legend()
          elif(degree == 5):
              five.append(test_pred)
              plt.scatter(100, test_pred2, color='turquoise')
              plt.plot(x_simulate, test_pred, color='turquoise', label = 'Degree 5')
              plt.legend()
          elif(degree == 6):
              six.append(test_pred)
              plt.scatter(100, test_pred2, color='brown')
              plt.plot(x_simulate, test_pred, color='brown', label = 'Degree 6')
              plt.legend()
          elif(degree == 7):
              seven.append(test_pred)
              plt.scatter(100, test_pred2, color='orange')
              plt.plot(x_simulate, test_pred, color='orange', label = 'Degree 7')
              plt.legend()
          else:
              eight.append(test_pred)
              plt.scatter(100, test_pred2, color='red')
              plt.plot(x_simulate, test_pred, color='red', label = 'Degree 8')
              plt.legend()'''

      regr_1 = DecisionTreeRegressor(max_depth=2)
      #regr_2 = DecisionTreeRegressor(max_depth=5)
      regr_1.fit(X, y_formatted)
      #regr_2.fit(X, y_formatted)

      # Predict
      test_predy_1 = regr_1.predict(x_simulate)
      #test_predy_2 = regr_2.predict(x_simulate)

      y_1 = regr_1.predict(X_prediction)
      #y_2 = regr_2.predict(X_prediction)

      # Plot the results
      plt.scatter(X, y_formatted, s=20, edgecolor="black",c="darkorange", label="data")
      plt.scatter(100,y_1,color="orange")
      #plt.scatter(100,y_2, color = "red")
      plt.plot(x_simulate, test_predy_1, color='orange', label = 'Depth = 4')
      #plt.plot(x_simulate, test_predy_2, color = 'red', label = 'Depth = 5')
      plt.legend()

   def percent_error(number):
       for i in range(len(finals)):
           number_indiviual = number[i-1]
           final_individual = finals[i-1]
           error_indiviual = np.abs(number_indiviual - final_individual) / final_individual * 100
           error.append(error_indiviual)
       calc.append(np.mean(error)/len(error))
       if(len(calc) == (17 or 16)):
           print(calc[len(calc) - 1])



   def log_regression(plt, X, Y, x_simulate, target):
      # create logarithmic regression model (didn't work)
      #logr = linear_model.LogisticRegression()
      #logr.fit(X, Y)
      #logyhat = logr.predict(X)
      #plt.plot(X, logyhat, color='orange')
      x = np.array(X, dtype=float) # transform your data in a numpy array of floats
      y = np.array(Y, dtype=float) # so the curve_fit can work
      try:
         popt, pcov = curve_fit(func_log, x.ravel(), y.ravel())
         #print ("a = %s , b = %s, c = %s" % (popt[0], popt[1], popt[2]))
         # print("Standard errors: ", np.sqrt(np.diag(pcov)))
         # predict value using logistic equation
         #target_memory = func_log(target, popt[0], popt[1], popt[2], popt[3])
         target_memory = func_log(target, *popt)
         plt.scatter(target, target_memory, color='purple')
         #plt.plot(x, func_log(x, *popt), color='purple')
         #y_simulate_log = func_log(x_simulate, popt[0], popt[1], popt[2], popt[3])
         y_simulate_log = func_log(x_simulate, *popt)
         plt.plot(x_simulate, y_simulate_log, color='purple', label = 'Logarithmic')
         plt.legend()
         predict_log = func_log(x,*popt)
         #print(r2_score(y,predict_log))
      except:
          popt, pcov = curve_fit(func_log2, x.ravel(), y.ravel())
          #print ("a = %s , b = %s, c = %s" % (popt[0], popt[1], popt[2]))
          # print("Standard errors: ", np.sqrt(np.diag(pcov)))
          # predict value using logistic equation
          #target_memory = func_log(target, popt[0], popt[1], popt[2], popt[3])
          target_memory = func_log2(target, *popt)
          plt.scatter(target, target_memory, color='purple')
          #plt.plot(x, func_log(x, *popt), color='purple')
          #y_simulate_log = func_log(x_simulate, popt[0], popt[1], popt[2], popt[3])
          y_simulate_log = func_log2(x_simulate, *popt)
          plt.plot(x_simulate, y_simulate_log, color='purple', label = 'Logarithmic')
          plt.legend()
          predict_log = func_log2(x,*popt)
          #print(r2_score(y,predict_log))
          #percent_error = (target_memory - final_memory)/final_memory * 100
          #print(task, "Percent Error : " + percent_error)
         # ignore RuntimeError: Optimal parameters not found: Number of calls to function has reached maxfev = 1000
   # log modeling
   def func_log(x, a, b, c, d):
      #return a * np.exp(-b * x) + c
      #return a * np.log(b * x) + c
      return a * np.log(b * (x + c)) + d
   def func_log2(x, a, b, d):
      #return a * np.exp(-b * x) + c
      #return a * np.log(b * x) + c
      return a * np.log(b * (x)) + d

   def get_memory_usage(result_dir, workflow_id, task, bucket):
      """ download task monitoring.log and get the max memory usage from the last line
         Args:
               result_dir: google cloud bucket workflow result directory
               workflow_id: workflow execution id
               task: task name
               bucket: google bucket
         Returns: a float number for task memory usage
      """
      base = result_dir + workflow_id + '/call-' + task
      try:
         path = base + '/shard-0/monitoring.log'
         blob = bucket.get_blob(path)
         blobBytes = blob.download_as_string()
      except AttributeError: # in case there's no shard folders
         path = base + '/monitoring.log'
         blob = bucket.get_blob(path)
         blobBytes = blob.download_as_string()
      blobString = blobBytes.decode("utf-8")
      index1 = find_second_last(blobString, '\n') + 1
      index2 = len(blobString) - 1
      return float(blobString[index1:index2])
   def find_second_last(text, pattern):
      return text.rfind(pattern, 0, text.rfind(pattern))

   def regression(known_data, target_value, filename='plot/hummingbird.png'):
      """Use sklearn regression module to predict target value and plot.
         Args:
               known_data: profiling data, used for regression
               target_value: value for target key of Downsample in configuration file
               filename: name of file to output
         Returns: an array of maximum values for each profile
      """
      X = np.array(list(known_data.keys())) # list is required for python 3
      # puts all the keys into a 2d array of [len(known_data)][1]
      X = X.reshape((len(known_data), 1))
      ys = np.array(list(known_data.values())) # list is required for python 3
      result = []
      x_test = np.linspace(1000, target_value, num=100)
      x_test = x_test.reshape(100, 1)
      plt.figure(figsize=(18, 10))
      for i in range(len(ys[0])):
         y = ys[:,i]
         regr = linear_model.LinearRegression()
         regr.fit(X, y)
         # Reshape data using array.reshape(-1, 1) if your data has a single feature
         target = np.array(target_value).reshape(-1, 1)
         # get maximum value of the maximum y-value vs the predicted y-value
         res = max(np.max(y), np.asscalar(regr.predict(target)))
         result.append(res)
         # plot data onto graph
         plt.subplot(len(ys[0]), 1, i + 1)
         plt.plot(x_test, regr.predict(x_test), 'b', X, y, 'ro')
      # save graph as a file
      plt.savefig(filename)
      return result
   if __name__== "__main__":
      main()
