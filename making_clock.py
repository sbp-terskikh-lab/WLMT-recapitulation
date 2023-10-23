# In[1]:

import pandas as pd
from time import gmtime, strftime
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
matplotlib.style.use('ggplot')
get_ipython().magic('matplotlib inline')
get_ipython().magic('pylab inline')
from sklearn import model_selection, linear_model, metrics
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline


# In[2]:

import seaborn
import sklearn

def plot_scores(optimizer):
    scores = [[item[0]['alpha'], 
               item[1], 
               (np.sum((item[2]-item[1])**2)/(item[2].size-1))**0.5] for item in optimizer.grid_scores_]
    scores = np.array(scores)
    plt.semilogx(scores[:,0], scores[:,1])
    plt.fill_between(scores[:,0], scores[:,1]-scores[:,2], 
                                  scores[:,1]+scores[:,2], alpha=0.3)
    plt.show()
    
def plot_scores_pipeline(optimizer):
    scores = [[item[0]['regression__alpha'], 
               item[1], 
               (np.sum((item[2]-item[1])**2)/(item[2].size-1))**0.5] for item in optimizer.grid_scores_]
    scores = np.array(scores)
    plt.semilogx(scores[:,0], scores[:,1])
    plt.fill_between(scores[:,0], scores[:,1]-scores[:,2], 
                                  scores[:,1]+scores[:,2], alpha=0.3)
    plt.show()
    
def heatmap_scores(optimizer,parameters_grid_alpha,parameters_grid_l1_ratio):
    
    scores = [[item[0]['alpha'],
           item[0]['l1_ratio'],
           item[1], 
           (np.sum((item[2]-item[1])**2)/(item[2].size-1))**0.5] for item in optimizer.grid_scores_]
    scores = np.array(scores)
    temp = pd.DataFrame(index=parameters_grid_l1_ratio,
                    columns=parameters_grid_alpha)
    for x in range(scores.shape[0]):
        temp.loc[scores[x,1],scores[x,0]] = scores[x,2]
    temp = temp.fillna(min(temp.min(0)))
    plt.figure(figsize=(8,6))
    seaborn.set(font_scale=1.5)
    cg = seaborn.heatmap(temp,cmap="Blues",linewidths=0.15)
    plt.yticks(rotation=0)
    cg.set(xlabel="alpha",ylabel="l1 ratio")
    
def heatmap_scores_pipeline(optimizer,parameters_grid_alpha,parameters_grid_l1_ratio):
    
    scores = [[item[0]['regression__alpha'],
           item[0]['regression__l1_ratio'],
           item[1], 
           (np.sum((item[2]-item[1])**2)/(item[2].size-1))**0.5] for item in optimizer.grid_scores_]
    scores = np.array(scores)
    temp = pd.DataFrame(index=parameters_grid_l1_ratio,
                    columns=parameters_grid_alpha)
    for x in range(scores.shape[0]):
        temp.loc[scores[x,1],scores[x,0]] = scores[x,2]
    temp = temp.fillna(min(temp.min(0)))
    plt.figure(figsize=(8,6))
    seaborn.set(font_scale=1.5)
    cg = seaborn.heatmap(temp,cmap="Blues",linewidths=0.15)
    plt.yticks(rotation=0)
    cg.set(xlabel="alpha",ylabel="l1 ratio")
    
def to_fraction(row):
    return(round(row/100,6))

def round_to(row,i):
    return(round(row,i))

def get_type(row):
    return(type(row))


# ____

# # Making clock

# In[37]:

current_random_state = 23
current_n_folds = 10
current_regression_alpha = [100]
current_regression_l1_ratio = [0.6]

################################
from sklearn.model_selection import train_test_split

(X_train, 
 X_test, 
 y_train, y_test) = train_test_split(metlev, ageslist, 
                                     test_size=0.2, 
                                     random_state=current_random_state,stratify=tissueslist)
################################

cv = model_selection.StratifiedKFold(y_train, n_folds = current_n_folds, shuffle = False, random_state = 23)
################################

#scaler = StandardScaler()
pensf_regressor = linear_model.ElasticNet(random_state=23)
pensf_pipeline = Pipeline(steps = [('regression', pensf_regressor)])
#pensf_pipeline.get_params()
scorer = metrics.make_scorer(metrics.mean_absolute_error, greater_is_better = False)
################################

parameters_grid = {
    'regression__alpha' : current_regression_alpha,
    #'regression__epsilon': [0.01,0.1,0.5],
    'regression__l1_ratio': current_regression_l1_ratio,
}  
################################

grid_pensf = model_selection.GridSearchCV(pensf_pipeline,parameters_grid,cv=cv,scoring=scorer,n_jobs=10)
################################



# In[38]:

grid_pensf.fit(X_train,y_train)
################################



# In[39]:

print(np.sum(grid_pensf.best_estimator_.named_steps['regression'].coef_!=0))
print(metrics.mean_absolute_error(y_test, grid_pensf.predict(X_test)))
print(metrics.mean_absolute_error(y_train, grid_pensf.predict(X_train)))
print(grid_pensf.best_params_)


# In[ ]:

y_true = y_train
y_pred = grid_pensf.best_estimator_.predict(X_train)
print(sklearn.metrics.r2_score(y_true, y_pred))


# In[ ]:

y_true = y_test
y_pred = grid_pensf.best_estimator_.predict(X_test)
print(sklearn.metrics.r2_score(y_true, y_pred))


# # Getting positions and weights

# In[77]:

potential_sites = grid_pensf.best_estimator_.named_steps['regression'].coef_
pos_number = []
weight = []
for i in range(len(potential_sites)):
    if potential_sites[i] != 0:
        pos_number.append(i)
        weight.append(potential_sites[i])

clock_positions = []
for i in range(len(wid_list)):
    if i in pos_number:
        clock_positions.append(wid_list[i])
        


# In[82]:

clock435 = pd.DataFrame()
clock435['index'] = clock_positions
clock435['weight'] = weight
clock435.head(10)


# In[ ]:

clock435.to_csv(WLMT_clock)

