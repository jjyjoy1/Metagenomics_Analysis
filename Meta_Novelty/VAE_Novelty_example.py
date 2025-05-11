# 1. Load k-mer/gene abundance matrix
import pandas as pd
data = pd.read_csv("kmer_abundance_matrix.csv", index_col=0)

# 2. Normalize and split
from sklearn.preprocessing import StandardScaler
X = StandardScaler().fit_transform(data)

# 3. Build VAE
from tensorflow.keras import layers, Model

latent_dim = 10
inputs = layers.Input(shape=(X.shape[1],))
h = layers.Dense(128, activation='relu')(inputs)
z_mean = layers.Dense(latent_dim)(h)
z_log_var = layers.Dense(latent_dim)(h)

# Sampling function
def sampling(args):
    z_mean, z_log_var = args
    epsilon = tf.random.normal(shape=tf.shape(z_mean))
    return z_mean + tf.exp(0.5 * z_log_var) * epsilon

z = layers.Lambda(sampling)([z_mean, z_log_var])

# Decoder
decoder_h = layers.Dense(128, activation='relu')(z)
outputs = layers.Dense(X.shape[1], activation='sigmoid')(decoder_h)

vae = Model(inputs, outputs)
vae.compile(optimizer='adam', loss='mse')
vae.fit(X, X, epochs=100, batch_size=32, validation_split=0.1)

# 4. Novelty scoring
X_pred = vae.predict(X)
recon_error = ((X - X_pred)**2).mean(axis=1)

# 5. Visualization
import matplotlib.pyplot as plt
import seaborn as sns
sns.histplot(recon_error, bins=50)

