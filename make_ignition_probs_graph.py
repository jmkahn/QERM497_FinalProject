#%%
import pandas as pd 
import matplotlib.pyplot as plt
ig_prob_df = pd.read_csv("outputs\ignition_probability_final.csv")

#%%
ignition_probability = ig_prob_df['ignition_probability']
min_biomass = ig_prob_df['min_biomass']
max_biomass = ig_prob_df['max_biomass']
mean_biomass = ig_prob_df['mean_biomass']

# Biomass: min, max, and mean. 
plt.figure()
plt.scatter(ignition_probability, max_biomass, color = "green", label = "max")
plt.scatter(ignition_probability, mean_biomass, color = "orange", label = "mean")
plt.scatter(ignition_probability, min_biomass, color = "blue", label = "min")
plt.legend()
plt.xlabel("ignition probability")
plt.ylabel("mean biomass")
plt.savefig("outputs/ignition_prob_biomass.png")
