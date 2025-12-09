"""
QUICK_START.py - Gene Regulation Analysis Tool
================================================
This script performs comprehensive gene regulation analysis including:
- Loading RNA-seq data from Excel files
- Identifying upregulated and downregulated genes
- Generating statistical analysis and visualizations
- Creating publication-ready figures

Author: RNAsq Team
Date: 2025-12-09
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Configure visualization style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)
plt.rcParams['font.size'] = 11


class GeneRegulationAnalyzer:
    """
    A comprehensive class for analyzing gene regulation patterns
    from RNA-seq data.
    """
    
    def __init__(self, excel_file, log2fc_col='log2FoldChange', 
                 padj_col='padj', gene_col='gene_name'):
        """
        Initialize the analyzer with RNA-seq data.
        
        Parameters:
        -----------
        excel_file : str
            Path to Excel file containing RNA-seq results
        log2fc_col : str
            Column name for log2 fold change values
        padj_col : str
            Column name for adjusted p-values
        gene_col : str
            Column name for gene names/identifiers
        """
        self.excel_file = excel_file
        self.log2fc_col = log2fc_col
        self.padj_col = padj_col
        self.gene_col = gene_col
        self.df = None
        self.upregulated = None
        self.downregulated = None
        self.stats_summary = {}
        
        print(f"Loading data from: {excel_file}")
        self._load_data()
    
    def _load_data(self):
        """Load and validate RNA-seq data from Excel file."""
        try:
            self.df = pd.read_excel(self.excel_file)
            print(f"✓ Data loaded successfully")
            print(f"  Total genes: {len(self.df)}")
            print(f"  Columns: {', '.join(self.df.columns.tolist())}")
        except FileNotFoundError:
            print(f"✗ Error: File '{self.excel_file}' not found")
            raise
        except Exception as e:
            print(f"✗ Error loading file: {e}")
            raise
    
    def identify_regulated_genes(self, log2fc_threshold=1.0, padj_threshold=0.05):
        """
        Identify upregulated and downregulated genes based on thresholds.
        
        Parameters:
        -----------
        log2fc_threshold : float
            Absolute log2 fold change threshold (default: 1.0 = 2-fold)
        padj_threshold : float
            Adjusted p-value significance threshold (default: 0.05)
        
        Returns:
        --------
        dict : Summary statistics of regulation
        """
        # Apply filtering criteria
        significant = (self.df[self.padj_col] < padj_threshold)
        
        self.upregulated = self.df[
            (self.df[self.log2fc_col] > log2fc_threshold) & significant
        ].copy()
        
        self.downregulated = self.df[
            (self.df[self.log2fc_col] < -log2fc_threshold) & significant
        ].copy()
        
        # Calculate statistics
        self.stats_summary = {
            'total_genes': len(self.df),
            'significant_genes': significant.sum(),
            'upregulated': len(self.upregulated),
            'downregulated': len(self.downregulated),
            'log2fc_threshold': log2fc_threshold,
            'padj_threshold': padj_threshold,
            'upregulated_fold_change_mean': self.upregulated[self.log2fc_col].mean(),
            'downregulated_fold_change_mean': self.downregulated[self.log2fc_col].mean(),
        }
        
        return self._print_summary()
    
    def _print_summary(self):
        """Print summary statistics."""
        print("\n" + "="*60)
        print("GENE REGULATION ANALYSIS SUMMARY")
        print("="*60)
        print(f"Total genes analyzed: {self.stats_summary['total_genes']:,}")
        print(f"Significant genes (padj < {self.stats_summary['padj_threshold']}): {self.stats_summary['significant_genes']:,}")
        print(f"\n📈 UPREGULATED GENES: {self.stats_summary['upregulated']:,}")
        print(f"   Mean log2FC: {self.stats_summary['upregulated_fold_change_mean']:.2f}")
        print(f"\n📉 DOWNREGULATED GENES: {self.stats_summary['downregulated']:,}")
        print(f"   Mean log2FC: {self.stats_summary['downregulated_fold_change_mean']:.2f}")
        print("="*60 + "\n")
        
        return self.stats_summary
    
    def plot_volcano(self, output_file='volcano_plot.png', 
                     log2fc_threshold=1.0, padj_threshold=0.05):
        """
        Create a volcano plot showing fold change vs p-value.
        
        Parameters:
        -----------
        output_file : str
            Path to save the figure
        log2fc_threshold : float
            Fold change threshold for visualization
        padj_threshold : float
            P-value threshold for visualization
        """
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Create color categories
        colors = ['gray'] * len(self.df)
        for i, row in self.df.iterrows():
            if row[self.padj_col] < padj_threshold:
                if row[self.log2fc_col] > log2fc_threshold:
                    colors[i] = '#e74c3c'  # Red for upregulated
                elif row[self.log2fc_col] < -log2fc_threshold:
                    colors[i] = '#3498db'  # Blue for downregulated
        
        # Create scatter plot
        scatter = ax.scatter(
            self.df[self.log2fc_col],
            -np.log10(self.df[self.padj_col]),
            c=colors, alpha=0.6, s=30, edgecolors='black', linewidth=0.5
        )
        
        # Add threshold lines
        ax.axvline(x=log2fc_threshold, color='red', linestyle='--', 
                   linewidth=2, alpha=0.5, label=f'log2FC = ±{log2fc_threshold}')
        ax.axvline(x=-log2fc_threshold, color='red', linestyle='--', 
                   linewidth=2, alpha=0.5)
        ax.axhline(y=-np.log10(padj_threshold), color='green', linestyle='--', 
                   linewidth=2, alpha=0.5, label=f'padj = {padj_threshold}')
        
        # Labels and formatting
        ax.set_xlabel('log2 Fold Change', fontsize=12, fontweight='bold')
        ax.set_ylabel('-log10(adjusted p-value)', fontsize=12, fontweight='bold')
        ax.set_title('Volcano Plot: Gene Expression Changes', 
                     fontsize=14, fontweight='bold', pad=20)
        ax.legend(fontsize=10, loc='upper right')
        ax.grid(True, alpha=0.3)
        
        # Add gene counts
        up_count = self.stats_summary.get('upregulated', 0)
        down_count = self.stats_summary.get('downregulated', 0)
        textstr = f'Upregulated: {up_count}\nDownregulated: {down_count}'
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=11,
                verticalalignment='top', bbox=dict(boxstyle='round', 
                facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Volcano plot saved to: {output_file}")
        plt.close()
    
    def plot_ma(self, output_file='ma_plot.png', log2fc_threshold=1.0, padj_threshold=0.05):
        """
        Create an MA plot (log fold change vs average expression).
        
        Parameters:
        -----------
        output_file : str
            Path to save the figure
        log2fc_threshold : float
            Fold change threshold
        padj_threshold : float
            P-value threshold
        """
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Calculate M (log fold change) and A (average log expression)
        # Assuming there are count columns - adjust as needed
        M = self.df[self.log2fc_col]
        A = np.log10(self.df.iloc[:, 2:].mean(axis=1) + 1)  # Approximate A value
        
        # Color assignment
        colors = ['gray'] * len(self.df)
        for i, row in self.df.iterrows():
            if row[self.padj_col] < padj_threshold:
                if row[self.log2fc_col] > log2fc_threshold:
                    colors[i] = '#e74c3c'
                elif row[self.log2fc_col] < -log2fc_threshold:
                    colors[i] = '#3498db'
        
        # Create scatter plot
        ax.scatter(A, M, c=colors, alpha=0.6, s=30, edgecolors='black', linewidth=0.5)
        
        # Add threshold lines
        ax.axhline(y=log2fc_threshold, color='red', linestyle='--', 
                   linewidth=2, alpha=0.5, label=f'log2FC = ±{log2fc_threshold}')
        ax.axhline(y=-log2fc_threshold, color='red', linestyle='--', 
                   linewidth=2, alpha=0.5)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.3)
        
        ax.set_xlabel('A = log10(Average Expression)', fontsize=12, fontweight='bold')
        ax.set_ylabel('M = log2(Fold Change)', fontsize=12, fontweight='bold')
        ax.set_title('MA Plot: Expression Magnitude vs Change', 
                     fontsize=14, fontweight='bold', pad=20)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ MA plot saved to: {output_file}")
        plt.close()
    
    def plot_distribution(self, output_file='distribution_plot.png'):
        """
        Create distribution plots for fold changes and p-values.
        
        Parameters:
        -----------
        output_file : str
            Path to save the figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # log2FC distribution
        axes[0, 0].hist(self.df[self.log2fc_col], bins=50, color='skyblue', 
                        edgecolor='black', alpha=0.7)
        axes[0, 0].axvline(x=0, color='red', linestyle='--', linewidth=2)
        axes[0, 0].set_xlabel('log2 Fold Change', fontweight='bold')
        axes[0, 0].set_ylabel('Frequency', fontweight='bold')
        axes[0, 0].set_title('Distribution of log2FC', fontweight='bold')
        axes[0, 0].grid(True, alpha=0.3)
        
        # log10(p-value) distribution
        axes[0, 1].hist(-np.log10(self.df[self.padj_col]), bins=50, 
                        color='lightcoral', edgecolor='black', alpha=0.7)
        axes[0, 1].set_xlabel('-log10(adjusted p-value)', fontweight='bold')
        axes[0, 1].set_ylabel('Frequency', fontweight='bold')
        axes[0, 1].set_title('Distribution of -log10(padj)', fontweight='bold')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Fold change box plot
        data_for_box = [
            self.df[self.log2fc_col],
            self.upregulated[self.log2fc_col] if len(self.upregulated) > 0 else [],
            self.downregulated[self.log2fc_col] if len(self.downregulated) > 0 else []
        ]
        bp = axes[1, 0].boxplot(data_for_box, labels=['All Genes', 'Upregulated', 'Downregulated'],
                                patch_artist=True)
        colors_box = ['skyblue', '#e74c3c', '#3498db']
        for patch, color in zip(bp['boxes'], colors_box):
            patch.set_facecolor(color)
        axes[1, 0].set_ylabel('log2 Fold Change', fontweight='bold')
        axes[1, 0].set_title('Fold Change by Regulation Status', fontweight='bold')
        axes[1, 0].grid(True, alpha=0.3, axis='y')
        
        # Gene count bar plot
        categories = ['All Genes', 'Upregulated', 'Downregulated']
        counts = [
            len(self.df),
            len(self.upregulated),
            len(self.downregulated)
        ]
        bars = axes[1, 1].bar(categories, counts, color=['skyblue', '#e74c3c', '#3498db'],
                              edgecolor='black', alpha=0.7)
        axes[1, 1].set_ylabel('Number of Genes', fontweight='bold')
        axes[1, 1].set_title('Gene Count Summary', fontweight='bold')
        axes[1, 1].grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            axes[1, 1].text(bar.get_x() + bar.get_width()/2., height,
                           f'{int(height):,}', ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Distribution plots saved to: {output_file}")
        plt.close()
    
    def export_regulated_genes(self, output_prefix='regulated_genes'):
        """
        Export upregulated and downregulated genes to separate files.
        
        Parameters:
        -----------
        output_prefix : str
            Prefix for output file names
        """
        # Export upregulated genes
        upregulated_file = f'{output_prefix}_upregulated.csv'
        self.upregulated.to_csv(upregulated_file, index=False)
        print(f"✓ Upregulated genes exported to: {upregulated_file}")
        
        # Export downregulated genes
        downregulated_file = f'{output_prefix}_downregulated.csv'
        self.downregulated.to_csv(downregulated_file, index=False)
        print(f"✓ Downregulated genes exported to: {downregulated_file}")
        
        # Export summary
        summary_file = f'{output_prefix}_summary.txt'
        with open(summary_file, 'w') as f:
            f.write("GENE REGULATION ANALYSIS SUMMARY\n")
            f.write("=" * 60 + "\n")
            for key, value in self.stats_summary.items():
                if isinstance(value, float):
                    f.write(f"{key}: {value:.4f}\n")
                else:
                    f.write(f"{key}: {value}\n")
        print(f"✓ Summary exported to: {summary_file}")


def main():
    """
    Main function demonstrating the complete analysis workflow.
    """
    print("\n" + "="*60)
    print("GENE REGULATION ANALYSIS - QUICK START")
    print("="*60 + "\n")
    
    # Configuration
    EXCEL_FILE = 'rna_seq_results.xlsx'  # Update with your file path
    LOG2FC_THRESHOLD = 1.0  # 2-fold change
    PADJ_THRESHOLD = 0.05
    
    try:
        # Initialize analyzer
        analyzer = GeneRegulationAnalyzer(
            excel_file=EXCEL_FILE,
            log2fc_col='log2FoldChange',  # Adjust column names as needed
            padj_col='padj',
            gene_col='gene_name'
        )
        
        # Identify regulated genes
        print("\n🔍 Identifying regulated genes...")
        analyzer.identify_regulated_genes(
            log2fc_threshold=LOG2FC_THRESHOLD,
            padj_threshold=PADJ_THRESHOLD
        )
        
        # Generate visualizations
        print("\n📊 Generating visualizations...")
        analyzer.plot_volcano(
            output_file='volcano_plot.png',
            log2fc_threshold=LOG2FC_THRESHOLD,
            padj_threshold=PADJ_THRESHOLD
        )
        
        analyzer.plot_ma(
            output_file='ma_plot.png',
            log2fc_threshold=LOG2FC_THRESHOLD,
            padj_threshold=PADJ_THRESHOLD
        )
        
        analyzer.plot_distribution(
            output_file='distribution_plots.png'
        )
        
        # Export results
        print("\n💾 Exporting results...")
        analyzer.export_regulated_genes(output_prefix='regulated_genes')
        
        print("\n✅ Analysis complete!")
        print("\nGenerated files:")
        print("  - volcano_plot.png")
        print("  - ma_plot.png")
        print("  - distribution_plots.png")
        print("  - regulated_genes_upregulated.csv")
        print("  - regulated_genes_downregulated.csv")
        print("  - regulated_genes_summary.txt\n")
        
    except FileNotFoundError:
        print(f"\n⚠️  Please ensure '{EXCEL_FILE}' exists in the current directory")
        print("Update the EXCEL_FILE variable with the correct path and try again.")
    except Exception as e:
        print(f"\n❌ Error during analysis: {e}")
        raise


if __name__ == "__main__":
    main()
