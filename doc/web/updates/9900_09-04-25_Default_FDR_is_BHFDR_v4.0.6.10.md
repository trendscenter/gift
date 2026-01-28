# Default FDR is BHFDR v4.0.6.10<BR>
Removed all instances of false discovery rate (FDR) only and replaced it with either MAFDR (Storey's FDR method) or BHFDR (Benjamini-Hochberg FDR method), which now is defauld. BHFDR may be used for smaller datasets with fewer p-values. MAFDR may be used in extremely large datasets with at least above 100 p-values.
