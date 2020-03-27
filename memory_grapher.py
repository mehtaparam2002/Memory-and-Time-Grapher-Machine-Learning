import warnings
import math
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
from pandas import DataFrame
from scipy.optimize import curve_fit
from sklearn.linear_model import LogisticRegression
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
import seaborn as sns
from sklearn import metrics



with warnings.catch_warnings():
   warnings.simplefilter("ignore")
   from google.cloud import storage
   def main():
      tasks = ['ataqc', 'bam2ta', 'bowtie2', 'filter', 'idr_pr',
      'macs2', 'macs2_pr1', 'macs2_pr2', 'macs2_signal_track', 'overlap_pr',
      'qc_report', 'read_genome_tsv', 'reproducibility_idr', 'reproducibility_overlap', 'spr',
      'trim_adapter', 'xcor']

      workflow_ids = [
         'b8d50752-4641-4a22-9875-effab47858f9',
         '762b3d68-5a7f-46eb-a895-ca2698796aa5',
         'c13c57ea-673b-4b30-bf6d-c8058eebc0fd']

      sizes = [
         0.1,
         0.1,
         1,
         1,
         2,
         5,
         10,
         20
      ]
      client = storage.Client()
      bucket = client.get_bucket('my-awesome-bucket-pmehta12')
      result_dir = 'ENCSR356KRQ_subsampled_short/atac/'
      for task in tasks:
         # create data frame
         memoryUsages = []
         for i, _ in enumerate(sizes):
            memoryUsages.append(getMemoryUsage(result_dir, workflow_ids[i], task, bucket))
         data = {
            'Size': sizes,
            'Memory': memoryUsages
         }
         plotMemoryUsage(data, task)
   def plotMemoryUsage(data, task):
      """ plot task memory usage of downsamples and predict memory usage of real sample
         Args:
               data: a table of sizes and memory usage arraya
               task: task name
         Returns: none
      """
      df = DataFrame(data, columns=['Size', 'Memory'])
      X = df[['Size']]
      Y = df['Memory']

      # create linear regression model
      regr = linear_model.LinearRegression()
      regr.fit(X, Y)
      yhat = regr.predict(X)
      # predict target memory usage. input should be a 2D array
      target = 100
      target_memory = regr.predict([[ target ]])
      prediction = 'real sample: ' + str(int(target_memory[0]))  + ' MB'
      # get r-squared
      score = '{:0.5f}'.format(regr.score(X, Y))
      r_squared =  '$R^2$ = ' + str(score)
      # get mean squared error
      mse = mean_squared_error(yhat, Y)
      rmse = 'RMSE = ' + str(int(math.sqrt(mse))) + ' MB'
      # get equation
      #coefficient = regr.coef_
      #intercept = regr.intercept_
      #equation = 'Y = ' + str(int(intercept)) + ' + ' + str(coefficient[0]) + 'X'
      ax = plt.subplot(1,1,1)
      plt.scatter(X, Y, color='red')
      plt.plot(X, yhat, color='blue')

      #Log regression

      plt.plot(x, y, 'ro',label="Original Data")
      x = np.array(x, dtype=float) #transform your data in a numpy array of floats
      y = np.array(y, dtype=float) #so the curve_fit can work

      def func(x, a, b, c, d):
          return a*x**3 + b*x**2 +c*x + d

      popt, pcov = curve_fit(func, x, y)

      print "(popt[0], popt[1], popt[2], popt[3])


      #Exponentional Regression


      plt.xlabel('Down Sample Size (%)')
      plt.ylabel('Memory Usage (MB)')
      plt.title(task)
      plt.text(0.5, 0.9, r_squared + '\n' + rmse + '\n' + prediction, ha='center', va='center', transform = ax.transAxes)
      plt.savefig('plot/' + task + '.png')
      plt.clf()
   def getMemoryUsage(result_dir, workflow_id, task, bucket):
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

if __name__== "__main__":
    main()
