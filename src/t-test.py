import numpy as np
from scipy.stats import pearsonr

# Data arrays
model1_remodelling_burden = np.array([64, 46, 48, 35, 50, 35, 28, 55])
model1_mPAP = np.array([32, 30, 22, 30, 26, 25, 27, 36])
model1_PVR = np.array([6.37, 3.45, 5.12, 2.57, 3.54, 3.85, 2.36, 8.61])

model2_remodelling_burden = np.array([59, 52, 55, 31, 47, 33, 32, 56])
model2_mPAP = np.array([19, 22, 16, 28, 22, 23, 24, 21])
model2_PVR = np.array([2.77, 2.0, 3.66, 2.26, 2.65, 3.46, 1.89, 4.41])

adjusted_scenario_mPAP = np.array([20, 26, 16, 21, 23, 22, 22, 24])
adjusted_scenario_PVR = np.array([2.72, 1.88, 2.77, 2.97, 2.57, 3.91, 1.97, 3.72])

# Correlation analyses for the provided scenarios
correlation_model1_mPAP_new = pearsonr(model1_remodelling_burden, model1_mPAP)
correlation_model1_PVR_new = pearsonr(model1_remodelling_burden, model1_PVR)
correlation_model2_mPAP_new = pearsonr(model2_remodelling_burden, model2_mPAP)
correlation_model2_PVR_new = pearsonr(model2_remodelling_burden, model2_PVR)

# Correlation for adjusted scenarios using model 2 remodelling burden
correlation_adjusted_mPAP = pearsonr(model2_remodelling_burden, adjusted_scenario_mPAP)
correlation_adjusted_PVR = pearsonr(model2_remodelling_burden, adjusted_scenario_PVR)

# Extracting p-values
p_values = {
    'Model 1 mPAP': correlation_model1_mPAP_new[1],
    'Model 1 PVR': correlation_model1_PVR_new[1],
    'Model 2 mPAP': correlation_model2_mPAP_new[1],
    'Model 2 PVR': correlation_model2_PVR_new[1],
    'Adjusted Scenario mPAP': correlation_adjusted_mPAP[1],
    'Adjusted Scenario PVR': correlation_adjusted_PVR[1]
}

print(p_values)

