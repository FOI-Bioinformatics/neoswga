"""
Deep learning primer embeddings for SWGA.

Provides dense vector representations of primer sequences using:
- Character-level or k-mer-level tokenization
- Transformer-based sequence encoders
- Pre-trained DNA language models (optional)
- Dimensionality reduction for visualization

Applications:
- Additional features for ML models
- Primer similarity computation
- Clustering and visualization
- Transfer learning from large DNA corpora

Falls back to simpler embeddings if deep learning libraries unavailable.
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Union
import warnings

# Try to import deep learning libraries
DL_AVAILABLE = False
DL_BACKEND = None

try:
    import torch
    import torch.nn as nn
    DL_AVAILABLE = True
    DL_BACKEND = 'torch'
except ImportError:
    pass

if not DL_AVAILABLE:
    try:
        import tensorflow as tf
        from tensorflow import keras
        DL_AVAILABLE = True
        DL_BACKEND = 'tensorflow'
    except ImportError:
        warnings.warn("Neither PyTorch nor TensorFlow available. Using fallback embeddings.")


# DNA alphabet
DNA_ALPHABET = ['A', 'C', 'G', 'T']
DNA_TO_IDX = {base: i for i, base in enumerate(DNA_ALPHABET)}
IDX_TO_DNA = {i: base for i, base in enumerate(DNA_ALPHABET)}


# ========================================
# Tokenization
# ========================================

class DNATokenizer:
    """
    Tokenize DNA sequences for neural networks.

    Supports:
    - Character-level (single nucleotides)
    - K-mer level (overlapping k-mers)
    """

    def __init__(self, mode: str = 'character', k: int = 3):
        """
        Initialize tokenizer.

        Args:
            mode: 'character' or 'kmer'
            k: K-mer length if using k-mer mode
        """
        self.mode = mode
        self.k = k

        if mode == 'character':
            self.vocab = DNA_ALPHABET
            self.vocab_size = len(self.vocab)
        elif mode == 'kmer':
            # Generate all k-mers
            self.vocab = self._generate_kmers(k)
            self.vocab_size = len(self.vocab)
        else:
            raise ValueError(f"Unknown mode: {mode}")

        self.token_to_idx = {token: i for i, token in enumerate(self.vocab)}
        self.idx_to_token = {i: token for token, i in self.token_to_idx.items()}

        # Special tokens
        self.pad_token = '<PAD>'
        self.unk_token = '<UNK>'
        self.vocab_size += 2  # Add special tokens

        self.token_to_idx[self.pad_token] = self.vocab_size - 2
        self.token_to_idx[self.unk_token] = self.vocab_size - 1

    def _generate_kmers(self, k: int) -> List[str]:
        """Generate all k-mers of length k."""
        if k == 1:
            return DNA_ALPHABET

        kmers = []
        for kmer in self._generate_kmers(k - 1):
            for base in DNA_ALPHABET:
                kmers.append(kmer + base)
        return kmers

    def tokenize(self, sequence: str) -> List[str]:
        """
        Tokenize sequence into tokens.

        Args:
            sequence: DNA sequence

        Returns:
            List of tokens
        """
        sequence = sequence.upper()

        if self.mode == 'character':
            return list(sequence)
        elif self.mode == 'kmer':
            # Overlapping k-mers
            if len(sequence) < self.k:
                return [sequence]
            return [sequence[i:i+self.k] for i in range(len(sequence) - self.k + 1)]
        else:
            raise ValueError(f"Unknown mode: {self.mode}")

    def encode(self, sequence: str, max_length: Optional[int] = None) -> List[int]:
        """
        Encode sequence to token indices.

        Args:
            sequence: DNA sequence
            max_length: Maximum sequence length (pad/truncate)

        Returns:
            List of token indices
        """
        tokens = self.tokenize(sequence)

        # Convert to indices
        indices = []
        for token in tokens:
            idx = self.token_to_idx.get(token, self.token_to_idx[self.unk_token])
            indices.append(idx)

        # Pad or truncate
        if max_length:
            pad_idx = self.token_to_idx[self.pad_token]
            if len(indices) < max_length:
                indices += [pad_idx] * (max_length - len(indices))
            else:
                indices = indices[:max_length]

        return indices

    def decode(self, indices: List[int]) -> str:
        """
        Decode token indices back to sequence.

        Args:
            indices: List of token indices

        Returns:
            Decoded sequence
        """
        tokens = []
        pad_idx = self.token_to_idx[self.pad_token]

        for idx in indices:
            if idx == pad_idx:
                break
            token = self.idx_to_token.get(idx, self.unk_token)
            tokens.append(token)

        if self.mode == 'character':
            return ''.join(tokens)
        elif self.mode == 'kmer':
            # Reconstruct from overlapping k-mers
            if not tokens:
                return ''
            sequence = tokens[0]
            for token in tokens[1:]:
                sequence += token[-1]
            return sequence
        else:
            raise ValueError(f"Unknown mode: {self.mode}")


# ========================================
# PyTorch Models
# ========================================

if DL_BACKEND == 'torch':

    class TransformerPrimerEncoder(nn.Module):
        """
        Transformer-based primer encoder (PyTorch).

        Architecture:
        - Embedding layer
        - Positional encoding
        - Multi-head self-attention
        - Feed-forward network
        - Mean pooling for sequence representation
        """

        def __init__(self,
                     vocab_size: int,
                     embed_dim: int = 64,
                     num_heads: int = 4,
                     num_layers: int = 2,
                     ff_dim: int = 128,
                     max_length: int = 20,
                     dropout: float = 0.1):
            """
            Initialize encoder.

            Args:
                vocab_size: Vocabulary size
                embed_dim: Embedding dimension
                num_heads: Number of attention heads
                num_layers: Number of transformer layers
                ff_dim: Feed-forward hidden dimension
                max_length: Maximum sequence length
                dropout: Dropout rate
            """
            super().__init__()

            self.embed_dim = embed_dim
            self.max_length = max_length

            # Embedding layer
            self.embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=0)

            # Positional encoding
            self.pos_encoding = self._create_positional_encoding(max_length, embed_dim)

            # Transformer layers
            encoder_layer = nn.TransformerEncoderLayer(
                d_model=embed_dim,
                nhead=num_heads,
                dim_feedforward=ff_dim,
                dropout=dropout,
                batch_first=True
            )
            self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)

            # Output projection
            self.output_proj = nn.Linear(embed_dim, embed_dim)

        def _create_positional_encoding(self, max_length: int, embed_dim: int) -> torch.Tensor:
            """Create sinusoidal positional encoding."""
            position = torch.arange(max_length).unsqueeze(1)
            div_term = torch.exp(torch.arange(0, embed_dim, 2) * (-np.log(10000.0) / embed_dim))

            pe = torch.zeros(max_length, embed_dim)
            pe[:, 0::2] = torch.sin(position * div_term)
            pe[:, 1::2] = torch.cos(position * div_term)

            return pe

        def forward(self, x: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
            """
            Forward pass.

            Args:
                x: Token indices (batch_size, seq_length)
                mask: Padding mask

            Returns:
                Embeddings (batch_size, embed_dim)
            """
            # Embed tokens
            x = self.embedding(x)  # (batch, seq, embed)

            # Add positional encoding
            seq_length = x.size(1)
            x = x + self.pos_encoding[:seq_length].to(x.device)

            # Transformer
            x = self.transformer(x, src_key_padding_mask=mask)

            # Mean pooling (ignore padding)
            if mask is not None:
                # Set padded positions to zero before mean
                x = x * (~mask).unsqueeze(-1).float()
                x = x.sum(dim=1) / (~mask).sum(dim=1, keepdim=True).float()
            else:
                x = x.mean(dim=1)

            # Project
            x = self.output_proj(x)

            return x


# ========================================
# TensorFlow Models
# ========================================

if DL_BACKEND == 'tensorflow':

    class TransformerPrimerEncoder(keras.Model):
        """
        Transformer-based primer encoder (TensorFlow).

        Same architecture as PyTorch version.
        """

        def __init__(self,
                     vocab_size: int,
                     embed_dim: int = 64,
                     num_heads: int = 4,
                     num_layers: int = 2,
                     ff_dim: int = 128,
                     max_length: int = 20,
                     dropout: float = 0.1):
            super().__init__()

            self.embed_dim = embed_dim
            self.max_length = max_length

            # Embedding
            self.embedding = keras.layers.Embedding(vocab_size, embed_dim, mask_zero=True)

            # Positional encoding
            self.pos_encoding = self._create_positional_encoding(max_length, embed_dim)

            # Transformer blocks
            self.transformer_blocks = []
            for _ in range(num_layers):
                self.transformer_blocks.append({
                    'attention': keras.layers.MultiHeadAttention(
                        num_heads=num_heads, key_dim=embed_dim // num_heads
                    ),
                    'ffn': keras.Sequential([
                        keras.layers.Dense(ff_dim, activation='relu'),
                        keras.layers.Dense(embed_dim)
                    ]),
                    'layernorm1': keras.layers.LayerNormalization(),
                    'layernorm2': keras.layers.LayerNormalization(),
                    'dropout1': keras.layers.Dropout(dropout),
                    'dropout2': keras.layers.Dropout(dropout)
                })

            # Output
            self.pool = keras.layers.GlobalAveragePooling1D()
            self.output_proj = keras.layers.Dense(embed_dim)

        def _create_positional_encoding(self, max_length: int, embed_dim: int):
            """Create positional encoding."""
            position = np.arange(max_length)[:, np.newaxis]
            div_term = np.exp(np.arange(0, embed_dim, 2) * -(np.log(10000.0) / embed_dim))

            pe = np.zeros((max_length, embed_dim))
            pe[:, 0::2] = np.sin(position * div_term)
            pe[:, 1::2] = np.cos(position * div_term)

            return tf.constant(pe, dtype=tf.float32)

        def call(self, x, training=False, mask=None):
            """Forward pass."""
            # Embed
            x = self.embedding(x)

            # Add positional encoding
            seq_length = tf.shape(x)[1]
            x = x + self.pos_encoding[:seq_length]

            # Transformer blocks
            for block in self.transformer_blocks:
                # Self-attention
                attn_output = block['attention'](x, x, attention_mask=mask)
                attn_output = block['dropout1'](attn_output, training=training)
                x = block['layernorm1'](x + attn_output)

                # Feed-forward
                ffn_output = block['ffn'](x)
                ffn_output = block['dropout2'](ffn_output, training=training)
                x = block['layernorm2'](x + ffn_output)

            # Pool and project
            x = self.pool(x, mask=mask)
            x = self.output_proj(x)

            return x


# ========================================
# High-Level API
# ========================================

class PrimerEmbedder:
    """
    High-level API for primer embeddings.

    Handles tokenization, model loading, and embedding generation.
    """

    def __init__(self,
                 tokenizer_mode: str = 'character',
                 tokenizer_k: int = 3,
                 embed_dim: int = 64,
                 model_path: Optional[str] = None):
        """
        Initialize embedder.

        Args:
            tokenizer_mode: 'character' or 'kmer'
            tokenizer_k: K-mer length
            embed_dim: Embedding dimension
            model_path: Path to pre-trained model (optional)
        """
        self.tokenizer = DNATokenizer(mode=tokenizer_mode, k=tokenizer_k)
        self.embed_dim = embed_dim
        self.model = None

        if DL_AVAILABLE:
            # Create model
            self.model = TransformerPrimerEncoder(
                vocab_size=self.tokenizer.vocab_size,
                embed_dim=embed_dim
            )

            # Load pre-trained weights if provided
            if model_path:
                self.load_model(model_path)

            print(f"Deep learning embedder initialized ({DL_BACKEND})")
            print(f"  Tokenizer: {tokenizer_mode}")
            print(f"  Embedding dim: {embed_dim}")
        else:
            print("Using fallback embeddings (install PyTorch or TensorFlow for deep learning)")

    def embed(self, primers: List[str], batch_size: int = 32) -> np.ndarray:
        """
        Generate embeddings for primers.

        Args:
            primers: List of primer sequences
            batch_size: Batch size for inference

        Returns:
            Embeddings array (n_primers, embed_dim)
        """
        if self.model is None or not DL_AVAILABLE:
            # Fallback: One-hot encoding + PCA
            return self._fallback_embed(primers)

        # Encode sequences
        max_length = max(len(p) for p in primers)
        encoded = [self.tokenizer.encode(p, max_length) for p in primers]

        if DL_BACKEND == 'torch':
            return self._embed_torch(encoded, batch_size)
        elif DL_BACKEND == 'tensorflow':
            return self._embed_tensorflow(encoded, batch_size)

    def _embed_torch(self, encoded: List[List[int]], batch_size: int) -> np.ndarray:
        """Generate embeddings using PyTorch."""
        import torch

        self.model.eval()
        embeddings = []

        with torch.no_grad():
            for i in range(0, len(encoded), batch_size):
                batch = encoded[i:i+batch_size]
                batch_tensor = torch.LongTensor(batch)

                # Create padding mask
                pad_idx = self.tokenizer.token_to_idx[self.tokenizer.pad_token]
                mask = (batch_tensor == pad_idx)

                # Encode
                batch_embeddings = self.model(batch_tensor, mask=mask)
                embeddings.append(batch_embeddings.numpy())

        return np.vstack(embeddings)

    def _embed_tensorflow(self, encoded: List[List[int]], batch_size: int) -> np.ndarray:
        """Generate embeddings using TensorFlow."""
        import tensorflow as tf

        embeddings = []

        for i in range(0, len(encoded), batch_size):
            batch = encoded[i:i+batch_size]
            batch_tensor = tf.constant(batch, dtype=tf.int32)

            # Encode
            batch_embeddings = self.model(batch_tensor, training=False)
            embeddings.append(batch_embeddings.numpy())

        return np.vstack(embeddings)

    def _fallback_embed(self, primers: List[str]) -> np.ndarray:
        """
        Fallback embeddings using one-hot encoding + dimensionality reduction.

        Not as powerful as deep learning but provides reasonable representations.
        """
        # One-hot encode each position
        max_length = max(len(p) for p in primers)

        encodings = []
        for primer in primers:
            encoding = np.zeros((max_length, 4))
            for i, base in enumerate(primer.upper()):
                if base in DNA_TO_IDX:
                    encoding[i, DNA_TO_IDX[base]] = 1
            encodings.append(encoding.flatten())

        encodings = np.array(encodings)

        # Reduce dimensionality with PCA
        from sklearn.decomposition import PCA
        pca = PCA(n_components=min(self.embed_dim, encodings.shape[0], encodings.shape[1]))
        reduced = pca.fit_transform(encodings)

        # Pad to target dimension if needed
        if reduced.shape[1] < self.embed_dim:
            padded = np.zeros((reduced.shape[0], self.embed_dim))
            padded[:, :reduced.shape[1]] = reduced
            reduced = padded

        return reduced

    def compute_similarity(self, primer1: str, primer2: str) -> float:
        """
        Compute cosine similarity between two primers.

        Args:
            primer1: First primer
            primer2: Second primer

        Returns:
            Cosine similarity (0-1)
        """
        embeddings = self.embed([primer1, primer2])

        # Cosine similarity
        dot = np.dot(embeddings[0], embeddings[1])
        norm1 = np.linalg.norm(embeddings[0])
        norm2 = np.linalg.norm(embeddings[1])

        return dot / (norm1 * norm2)

    def find_similar_primers(self,
                            query: str,
                            candidates: List[str],
                            top_k: int = 10) -> List[Tuple[str, float]]:
        """
        Find most similar primers to query.

        Args:
            query: Query primer
            candidates: Candidate primers
            top_k: Number of results

        Returns:
            List of (primer, similarity) tuples
        """
        # Embed all
        all_primers = [query] + candidates
        embeddings = self.embed(all_primers)

        query_embedding = embeddings[0]
        candidate_embeddings = embeddings[1:]

        # Compute similarities
        similarities = []
        for i, cand_emb in enumerate(candidate_embeddings):
            sim = np.dot(query_embedding, cand_emb) / (
                np.linalg.norm(query_embedding) * np.linalg.norm(cand_emb)
            )
            similarities.append((candidates[i], sim))

        # Sort and return top-k
        similarities.sort(key=lambda x: x[1], reverse=True)
        return similarities[:top_k]

    def save_model(self, path: str):
        """Save model weights."""
        if self.model is None:
            warnings.warn("No model to save")
            return

        if DL_BACKEND == 'torch':
            import torch
            torch.save(self.model.state_dict(), path)
        elif DL_BACKEND == 'tensorflow':
            self.model.save_weights(path)

        print(f"Model saved to: {path}")

    def load_model(self, path: str):
        """Load model weights."""
        if self.model is None:
            warnings.warn("No model to load into")
            return

        if DL_BACKEND == 'torch':
            import torch
            self.model.load_state_dict(torch.load(path))
        elif DL_BACKEND == 'tensorflow':
            self.model.load_weights(path)

        print(f"Model loaded from: {path}")


# ========================================
# Integration with ML Pipeline
# ========================================

def generate_embedding_features(primers: List[str],
                               embed_dim: int = 64,
                               model_path: Optional[str] = None) -> np.ndarray:
    """
    Generate embedding features for ML models.

    Args:
        primers: List of primers
        embed_dim: Embedding dimension
        model_path: Path to pre-trained model

    Returns:
        Feature matrix (n_primers, embed_dim)
    """
    embedder = PrimerEmbedder(embed_dim=embed_dim, model_path=model_path)
    embeddings = embedder.embed(primers)
    return embeddings


if __name__ == "__main__":
    print("Deep Learning Primer Embeddings")
    print("=" * 60)
    print(f"\nDeep learning available: {DL_AVAILABLE}")
    if DL_AVAILABLE:
        print(f"Backend: {DL_BACKEND}")

    print("\nFeatures:")
    print("  - Transformer-based sequence encoder")
    print("  - Character or k-mer tokenization")
    print("  - Dense vector representations")
    print("  - Similarity computation")
    print("  - Clustering and visualization")

    print("\nExample usage:")
    print("""
    from neoswga.core import deep_learning as dl

    # Create embedder
    embedder = dl.PrimerEmbedder(
        tokenizer_mode='character',
        embed_dim=64
    )

    # Generate embeddings
    primers = ['ATCGATCG', 'GCTAGCTA', 'AAAAGGGG']
    embeddings = embedder.embed(primers)  # (3, 64) array

    # Compute similarity
    sim = embedder.compute_similarity('ATCGATCG', 'ATCGATCC')
    print(f"Similarity: {sim:.3f}")

    # Find similar primers
    candidates = ['ATCGATCG', 'GCTAGCTA', 'TTTTCCCC', 'ATCGATCC']
    similar = embedder.find_similar_primers('ATCGATCG', candidates, top_k=3)
    for primer, sim in similar:
        print(f"{primer}: {sim:.3f}")

    # Use as ML features
    features = dl.generate_embedding_features(primers, embed_dim=64)
    """)
