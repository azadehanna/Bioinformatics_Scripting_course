#necessary libraries
import seaborn as sns
import numpy as np
import pandas as pd 
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.vq import kmeans2, whiten
from scipy.spatial.distance import cdist

#input files and directories
input_1= "./inputfiles/clinical_data.txt"
output="./clinical_data.stats.txt"
diversity_dir = "./inputfiles/diversityScores/"
distance_dir = "./inputfiles/distanceFiles/"

# Calculate Stat (Mean and Std)
def calculate_stat (animal):
    df=pd.read_table(animal)
    Mean_value=np.mean(df.iloc[:, 0])
    STD_value=np.std(df.iloc[:, 0])
    return  round(Mean_value,2), round(STD_value,2)

#add 2 stat values to clinical data

with open(input_1, "r") as file , open (output, "w") as output:
    header=file.readline().strip().split("\t")

    newheader= "\t".join(header+["Mean_value", "STD_value"])+"\n"
    output.write(newheader)

    for line in file:
        line_value=line.strip().split("\t")
        code_name=line_value[5]
        Mean_value , STD_value = calculate_stat (f"./inputfiles/diversityScores/{code_name}.diversity.txt")
        newline="\t".join(line_value+[str(Mean_value), str(STD_value)])+"\n"
        output.write(newline)


# read clinical data with stat
df_clinical = pd.read_csv("./clinical_data.stats.txt", sep="\t")
#print (df_clinical) 

# Sort to get top 2 and lowest 1 based on mean value
sorted_df = df_clinical.sort_values(by="Mean_value", ascending=False)
#print(sorted_df)
top_2 = sorted_df.iloc[:2]
lowest_1 = sorted_df.iloc[-1:]

#selected animals and print top2 and last one 
selected_animals = pd.concat([top_2, lowest_1])["code_name"].tolist()
print(f"Selected animals: {selected_animals}")    

def plot_and_cluster(animal, distance_dir):
    distance_file = os.path.join(distance_dir, f"{animal}.distance.txt")
    
    # Read the distance file
    df_animal = pd.read_csv(distance_file, header=None)
    
    # Create plot
    sns.set_style("whitegrid")  
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter plot with Seaborn
    sns.scatterplot(
        x=df_animal.iloc[:, 0], 
        y=df_animal.iloc[:, 1], 
        alpha=0.5, 
        ax=ax
    )
    
    # Set labels and title using the correct axis reference
    ax.set_title(f"Scatter Plot for {animal}", fontsize=14, fontweight='bold')
    ax.set_xlabel("X-axis Label", fontsize=12, fontweight='bold')
    ax.set_ylabel("Y-axis Label", fontsize=12, fontweight='bold')
    
    # Save the plot as both PNG and PDF
    plt.savefig(f"{animal}.png", format="png", bbox_inches="tight", pad_inches=0.5)
    plt.savefig(f"{animal}.pdf", format="pdf", bbox_inches="tight", pad_inches=0.5)
    plt.close()


# Loop through selected animals and generate plots
for animal in selected_animals:
    plot_and_cluster(animal, distance_dir)


#####Kmean 

for animal in selected_animals:
     distance_file = os.path.join(distance_dir, f"{animal}.distance.txt")
    
    # Read the distance file
     df_animal = pd.read_csv(distance_file, header=None)

# K-means clustering
     data = df_animal.to_numpy()
     data_whitened = whiten(data)

    # Determine optimal clusters using elbow method
     distortions = []
     K = range(1, 10)
     for k in K:
        centroids, labels = kmeans2(data_whitened, k, minit='points')
        distortions.append(sum(np.min(cdist(data_whitened, centroids, 'euclidean'), axis=1)) / data_whitened.shape[0])

    # Find elbow point
     optimal_k = 3  # Assuming elbow method determines 3 clusters
     plt.figure(figsize=(8, 6))
     plt.plot(K, distortions, 'bo-')
     plt.xlabel('Number of clusters')
     plt.ylabel('Distortion')
     plt.title(f'Elbow Method for {animal}')
     plt.savefig(f"{animal}_elbow.pdf")
     plt.close()

    # Perform K-means clustering
     centroids, labels = kmeans2(data_whitened, optimal_k, minit='points')

     plt.figure(figsize=(8, 6))
     sns.scatterplot(x=df_animal.iloc[:, 0], y=df_animal.iloc[:, 1], hue=labels, palette="viridis", alpha=0.6)
     plt.title(f"K-Means Clustering for {animal}")
     plt.xlabel("X-axis")
     plt.ylabel("Y-axis")
     plt.savefig(f"{animal}_clustered.pdf")
     plt.close()

