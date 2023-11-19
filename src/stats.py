import sklearn.metrics
from sklearn.metrics import confusion_matrix
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def calc_prec_recall(gt_labels, pred_labels, class_labels):
    '''
    Returns dictionary with prec:list, recall:list, f1:list, micro_prec:float, micro_recall:float, f1:float, acc:float
    '''
    if len(gt_labels) != len(pred_labels):
        raise ValueError("gt_labels and pred_labels must be same length")

    # ground truth are rows, predicted are columns
    #conf_matrix = confusion_matrix(gt_labels, pred_labels, labels = class_labels)
    p, r, fs, s = sklearn.metrics.precision_recall_fscore_support(gt_labels, pred_labels, labels=class_labels, zero_division=0.)

    stats = {}
    stats["prec"] = p
    stats["rec"] = r
    stats["f1"] = fs
    # these are all same for multiclass classification
    stats["micro_prec"] = sklearn.metrics.precision_score(gt_labels, pred_labels, labels=class_labels, average='micro', zero_division=0.)
    stats["micro_rec"] = sklearn.metrics.recall_score(gt_labels, pred_labels, labels=class_labels, average='micro', zero_division=0.)
    stats["micro_f1"] = sklearn.metrics.f1_score(gt_labels, pred_labels, labels=class_labels, average='micro', zero_division=0.)
    stats["acc"] = sklearn.metrics.accuracy_score(gt_labels, pred_labels, normalize=True)

    return stats

# i is class id:
# tp[i] = conf_matrix[i][i]
# tn[i] = conf_matrix_sum - rowsums[i] - colsums[i]
# fp[i] = colsums[i] - conf_matrix[i][i]
# fn[i] = rowsums[i] - conf_matrix[i][i]
# accuracy = diagsum / conf_matrix_sum
# precision[i] = tp[i]/(tp[i] + fp[i]) = np.trace(conf_matrix) / colsums
# recall[i] = tp[i]/(tp[i] + fn[i])    = np.trace(conf_matrix) / rowsums
# f1[i] = 2./( (1./precision[i]) + (1./recall[i])) = ...
# macro_precision = sum(precision) / (len(precision))
# macro_recall = sum(recall) / (len(recall))
#
# also, from https://www.evidentlyai.com/ml-in-production/data-drift:
# Data drift: is a change in the input data seen as a change in the data distribution
#    (presumably, the "bias" change as change in distribution location parameter, or "variance" change
#     as distribution scale parameter)
# Concept drift: is a change in input-output relationships, e.g. predictions are changing.
# model drift: models need to be maintained, that is retraining the models on the new data.
#      new data for example in products, customer preferences, and market condition changes.
#      - model monitoring: tracking metrics related to model quality (e.g., accuracy, precision, etc.),
#                           data and prediction drift, data quality, model bias, and fairness
#                https://www.evidentlyai.com/ml-in-production/model-monitoring
#         -- metrics for differenty models:
#            Classification: model accuracy, precision, recall, F1-score.
#            Regression: mean absolute error (MAE), mean squared error (MSE),
#                        mean absolute percentage error (MAPE), etc.
#            Ranking and recommendations: normalized discounted cumulative gain (NDCG),
#                        precision at K, mean average precision (MAP), etc.
# prediction drift:
# training-serving skew: mismatch between the data the model was trained on and the data it encounters in production.
# outliers

gt_labels_0=[]
gt_labels_0.extend([3,9])
gt_labels_0.extend([2, 2, 3, 9, 5, 2, 2, 2, 2, 2])
gt_labels_0.extend([3, 2, 8])
gt_labels_0.extend([8, 2, 8])
gt_labels_0.extend([ 8, 3, 8, 9, 9, 3])
gt_labels_0.extend([6])
gt_labels_0.extend([3, 5, 8, 9, 8, 0, 0, 9])
gt_labels_0.extend([4, 9, 9, 9, 9, 8, 0, 0])
gt_labels_0.extend([4, 0])
gt_labels_0.extend([4, 4, 9, 9, 0, 9, 0, 9, 0])
gt_labels_0.extend([3, 0, 3, 7])
gt_labels_0.extend([1, 6])
gt_labels_0.extend([9, 1])
gt_labels_0.extend([9, 6])

num_classes = 9 #ignoring dummy class id 9

scores_0= [0.29, 0.25, 0.55, 0.45, 0.39, 0.33, 0.32, 0.32, 0.30, 0.29, 0.25, 0.25, 0.53, 0.38, 0.29, 0.62, 0.32, 0.31, 0.66, 0.45, 0.43, 0.29, 0.28, 0.27, 0.32, 0.63, 0.47, 0.44, 0.42, 0.32, 0.29, 0.28, 0.26, 0.52, 0.48, 0.45, 0.37, 0.33, 0.30, 0.28, 0.28, 0.44, 0.34, 0.39, 0.36, 0.35, 0.33, 0.31, 0.31, 0.28, 0.27, 0.27, 0.61, 0.48, 0.35, 0.28, 0.35, 0.28, 0.50, 0.47, 0.38, 0.35]

pred_labels_0= [3, 7, 0, 0, 3, 7, 5, 1, 0, 2, 0, 2, 3, 1, 8, 8, 2, 1, 8, 3, 1, 4, 7, 5, 4, 3, 5, 8, 7, 1, 0, 0, 0, 7, 0, 0, 7, 0, 0, 7, 7, 7, 0, 8, 7, 5, 0, 0, 7, 0, 4, 7, 3, 7, 8, 2, 1, 7, 1, 1, 7, 6]

#dictionary with key = class id, value = dictionary with key=p or r, value = list
class_pr_lists = {} #dictionary key is class number
class_labels = []
for i in range(num_classes) :
    # each list element is from a threshold statistic
    class_pr_lists[i] = {'p':[], 'r':[]}
    class_labels.append(i)

micro_p = []
micro_r = []
micro_f1 = []
acc = []
_thresh = []

for thresh in np.arange(0.2, 0.95, 0.05):
    # indexes with scores >= thresh
    gt_labels = []
    pred_labels = []
    for i, s in enumerate(scores_0):
        if s >= thresh and (gt_labels_0[i] < 9):
            gt_labels.append(gt_labels_0[i])
            pred_labels.append(pred_labels_0[i])
    if len(gt_labels) == 0:
        break
    stats = calc_prec_recall(gt_labels, pred_labels, class_labels)

    p = stats["prec"]
    r = stats["rec"]
    for ii in range(num_classes) :
        class_pr_lists[ii]['p'].append(p[ii])
        class_pr_lists[ii]['r'].append(r[ii])
    micro_p.append(stats["micro_prec"])
    micro_r.append(stats["micro_rec"])
    micro_f1.append(stats["micro_f1"])
    acc.append(stats["acc"])
    _thresh.append(thresh)

ncols = 3
nrows = 3
fig = plt.figure(figsize=(6, 6))
#fig.suptitle("per class p vs r, thresh vs p,r,f1,acc")
# nrows, ncols, index where index starts at 1 in the upper left corner and increases to the right
#color='purple', marker='o', linestyle='solid', linewidth=2, markersize=1)

for i in range(num_classes):
    ax = plt.subplot(nrows, ncols, i+1)
    ax.set_title(str(i), x=0.2, y=0.8, color='red')
    plt.plot(class_pr_lists[i]['p'], class_pr_lists[i]['r'], 's', marker='o', markersize=2)
    if (i % 3) == 0:
        plt.ylabel('precision')
    if (i > 5):
        plt.xlabel("recall")
plt.show()

fig = plt.figure(figsize=(6, 6))
ax = plt.subplot(1, 1, 1)
plt.plot(_thresh, micro_p, label='p', marker='o', markersize=2, color='blue')
plt.plot(_thresh, micro_r, label='r', marker='o', markersize=2, color='red')
plt.plot(_thresh, micro_f1, label='f1', marker='o', markersize=2, color='orange')
plt.plot(_thresh, acc, label='a', marker='o', markersize=2, color='black', linestyle='dashed')
plt.xlabel("score threshhold")
plt.legend()
plt.show()
print(f'done')

'''
number in training set per class:
0 cupcake :              4
1 euclair:               4
2 icecream:              3
3 gingerbread_man:       5
4 icecream_sandwich:     3
5 honeycomb:             2
6 kitkat:                2
7 jellybean:             2
8 donut:                 2

for the 14 test frames, did the ground truth labels manually, lab_<img number>=[...].
used a label=9 for N/A as an identification not in classes (like a truck or an untrained statue)

lab_1=[3,9]
pred_1=[3, 7]
score_1= ['0.29', '0.25']

lab_2= [2, 2, 3, 9, 5, 2, 2, 2, 2, 2]
pred_2=[0, 0, 3, 7, 5, 1, 0, 2, 0, 2]
score_2= ['0.55', '0.45', '0.39', '0.33', '0.32', '0.32', '0.30', '0.29', '0.25', '0.25']

lab_3=[3, 2, 8]
pred_3=[3, 1, 8]
score_3= ['0.53', '0.38', '0.29']

lab_4=[8, 2, 8]
pred_4=[8, 2, 1]
score_4= ['0.62', '0.32', '0.31']

# possibly an error in labels here.  cannot see the title and don't want to re-run to print boxes
lab_5=[ 8, 3, 8, 9, 9, 3]
pred_5=[8, 3, 1, 4, 7, 5]
score_5= ['0.66', '0.45', '0.43', '0.29', '0.28', '0.27']

lab_6=[6]
pred_6=[4]
score_6= ['0.32']

lab_7 =[3, 5, 8, 9, 8, 0, 0, 9]
pred_7=[3, 5, 8, 7, 1, 0, 0, 0]
score_7= ['0.63', '0.47', '0.44', '0.42', '0.32', '0.29', '0.28', '0.26']

lab_8 =[4, 9, 9, 9, 9, 8, 0, 0]
pred_8=[7, 0, 0, 7, 0, 0, 7, 7]
score_8= ['0.52', '0.48', '0.45', '0.37', '0.33', '0.30', '0.28', '0.28']

lab_9=[4, 0]
pred_9=[7, 0]
pred_9=[7, 0]
score_9= ['0.44', '0.34']

lab_10 =[4, 4, 9, 9, 0, 9, 0, 9, 0]
pred_10=[8, 7, 5, 0, 0, 7, 0, 4, 7]
score_10= ['0.39', '0.36', '0.35', '0.33', '0.31', '0.31', '0.28', '0.27', '0.27']

lab_11=[3, 0, 3, 7]
pred_11=[3, 7, 8, 2]
score_11= ['0.61', '0.48', '0.35', '0.28']

lab_12=[1, 6]
pred_12=[1, 7]
score_12= ['0.35', '0.28']

lab_13=[9, 1]
pred_13=[1, 1]
score_13= ['0.50', '0.47']

lab_14=[9, 6]
pred_14=[7, 6]
score_14= ['0.38', '0.35']

pred_all= [3, 7, 0, 0, 3, 7, 5, 1, 0, 2, 0, 2, 3, 1, 8, 8, 2, 1, 8, 3, 1, 4, 7, 5, 4, 3, 5, 8, 7, 1, 0, 0, 0, 7, 0, 0, 7, 0, 0, 7, 7, 7, 0, 8, 7, 5, 0, 0, 7, 0, 4, 7, 3, 7, 8, 2, 1, 7, 1, 1, 7, 6]
scores_all= ['0.29', '0.25', '0.55', '0.45', '0.39', '0.33', '0.32', '0.32', '0.30', '0.29', '0.25', '0.25', '0.53', '0.38', '0.29', '0.62', '0.32', '0.31', '0.66', '0.45', '0.43', '0.29', '0.28', '0.27', '0.32', '0.63', '0.47', '0.44', '0.42', '0.32', '0.29', '0.28', '0.26', '0.52', '0.48', '0.45', '0.37', '0.33', '0.30', '0.28', '0.28', '0.44', '0.34', '0.39', '0.36', '0.35', '0.33', '0.31', '0.31', '0.28', '0.27', '0.27', '0.61', '0.48', '0.35', '0.28', '0.35', '0.28', '0.50', '0.47', '0.38', '0.35']

'''

