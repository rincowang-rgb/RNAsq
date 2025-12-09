#!/usr/bin/env python3
"""
Gene Regulation Analysis Script
================================

This script analyzes differential gene expression data to identify and characterize
upregulated and downregulated genes. It provides comprehensive statistical summaries
and generates publication-quality visualizations.

Features:
    - Loads gene expression data from Excel files
    - Separates upregulated from downregulated genes
    - Calculates summary statistics
    - Generates volcano plots and bar charts
    - Exports results to CSV files
    - Provides detailed logging and error handling

Author: RNAsq Analysis Pipeline
Date: 2025-12-09
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from typing import Tuple, Dict, List
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('gene_regulation_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class GeneRegulationAnalyzer:
    """
    A comprehensive analyzer for gene expression and regulation data.
    
    Attributes:
        data (pd.DataFrame): The loaded gene expression data
        input_file (str): Path to the input Excel file
        output_dir (Path): Directory for output files
        upregulated_genes (pd.DataFrame): Genes with positive Log2FoldChange
        downregulated_genes (pd.DataFrame): Genes with negative Log2FoldChange
    """
    
    def __init__(self, input_file: str, output_dir: str = 'results'):
        """
        Initialize the Gene Regulation Analyzer.
        
        Args:
            input_file (str): Path to the input Excel file
            output_dir (str): Directory to save output files (default: 'results')
        """
        self.input_file = input_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        self.data = None
        self.upregulated_genes = None
        self.downregulated_genes = None
        self.summary_stats = None
        
        logger.info(f"Gene Regulation Analyzer initialized with output directory: {self.output_dir}")
    
    def load_data(self) -> bool:
        """
        Load gene expression data from Excel file.
        
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            logger.info(f"Loading data from {self.input_file}...")
            self.data = pd.read_excel(self.input_file)
            logger.info(f"Successfully loaded {len(self.data)} genes")
            logger.info(f"Columns: {list(self.data.columns)}")
            return True
        except FileNotFoundError:
            logger.error(f"File not found: {self.input_file}")
            return False
        except Exception as e:
            logger.error(f"Error loading data: {str(e)}")
            return False
    
    def separate_genes(self, log2fc_column: str = 'Log2FoldChange', 
                      threshold: float = 0.0) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Separate genes into upregulated and downregulated based on Log2FoldChange.
        
        Args:
            log2fc_column (str): Name of the Log2FoldChange column (default: 'Log2FoldChange')
            threshold (float): Threshold for separation (default: 0.0)
        
        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: (upregulated_genes, downregulated_genes)
        """
        try:
            if log2fc_column not in self.data.columns:
                logger.error(f"Column '{log2fc_column}' not found in data")
                logger.info(f"Available columns: {list(self.data.columns)}")
                return None, None
            
            # Separate genes
            self.upregulated_genes = self.data[self.data[log2fc_column] > threshold].copy()
            self.downregulated_genes = self.data[self.data[log2fc_column] < threshold].copy()
            
            logger.info(f"Upregulated genes: {len(self.upregulated_genes)}")
            logger.info(f"Downregulated genes: {len(self.downregulated_genes)}")
            
            return self.upregulated_genes, self.downregulated_genes
        except Exception as e:
            logger.error(f"Error separating genes: {str(e)}")
            return None, None
    
    def display_top_genes(self, n: int = 10, log2fc_column: str = 'Log2FoldChange',
                         padj_column: str = 'padj') -> None:
        """
        Display and log the top upregulated and downregulated genes.
        
        Args:
            n (int): Number of top genes to display (default: 10)
            log2fc_column (str): Name of Log2FoldChange column (default: 'Log2FoldChange')
            padj_column (str): Name of adjusted p-value column (default: 'padj')
        """
        if self.upregulated_genes is None or self.downregulated_genes is None:
            logger.warning("Genes not separated yet. Call separate_genes() first.")
            return
        
        logger.info("=" * 80)
        logger.info(f"TOP {n} UPREGULATED GENES")
        logger.info("=" * 80)
        
        top_up = self.upregulated_genes.nlargest(n, log2fc_column)
        for idx, row in top_up.iterrows():
            gene_name = row.get('Gene', row.get('gene_id', 'Unknown'))
            log2fc = row.get(log2fc_column, np.nan)
            padj = row.get(padj_column, np.nan)
            logger.info(f"  {gene_name}: Log2FC={log2fc:.4f}, adj.p-value={padj:.2e}")
        
        logger.info("=" * 80)
        logger.info(f"TOP {n} DOWNREGULATED GENES")
        logger.info("=" * 80)
        
        top_down = self.downregulated_genes.nsmallest(n, log2fc_column)
        for idx, row in top_down.iterrows():
            gene_name = row.get('Gene', row.get('gene_id', 'Unknown'))
            log2fc = row.get(log2fc_column, np.nan)
            padj = row.get(padj_column, np.nan)
            logger.info(f"  {gene_name}: Log2FC={log2fc:.4f}, adj.p-value={padj:.2e}")
        
        logger.info("=" * 80)
    
    def calculate_summary_statistics(self, log2fc_column: str = 'Log2FoldChange',
                                    padj_column: str = 'padj') -> Dict:
        """
        Calculate comprehensive summary statistics.
        
        Args:
            log2fc_column (str): Name of Log2FoldChange column (default: 'Log2FoldChange')
            padj_column (str): Name of adjusted p-value column (default: 'padj')
        
        Returns:
            Dict: Dictionary containing summary statistics
        """
        if self.data is None or self.upregulated_genes is None or self.downregulated_genes is None:
            logger.warning("Data not properly loaded or separated")
            return None
        
        try:
            stats = {
                'total_genes': len(self.data),
                'upregulated_count': len(self.upregulated_genes),
                'downregulated_count': len(self.downregulated_genes),
                'upregulated_percent': (len(self.upregulated_genes) / len(self.data)) * 100,
                'downregulated_percent': (len(self.downregulated_genes) / len(self.data)) * 100,
                
                # Log2FoldChange statistics
                'upregulated_log2fc_mean': self.upregulated_genes[log2fc_column].mean(),
                'upregulated_log2fc_median': self.upregulated_genes[log2fc_column].median(),
                'upregulated_log2fc_std': self.upregulated_genes[log2fc_column].std(),
                'upregulated_log2fc_min': self.upregulated_genes[log2fc_column].min(),
                'upregulated_log2fc_max': self.upregulated_genes[log2fc_column].max(),
                
                'downregulated_log2fc_mean': self.downregulated_genes[log2fc_column].mean(),
                'downregulated_log2fc_median': self.downregulated_genes[log2fc_column].median(),
                'downregulated_log2fc_std': self.downregulated_genes[log2fc_column].std(),
                'downregulated_log2fc_min': self.downregulated_genes[log2fc_column].min(),
                'downregulated_log2fc_max': self.downregulated_genes[log2fc_column].max(),
            }
            
            # Add p-value statistics if available
            if padj_column in self.data.columns:
                stats.update({
                    'upregulated_padj_mean': self.upregulated_genes[padj_column].mean(),
                    'upregulated_padj_median': self.upregulated_genes[padj_column].median(),
                    'downregulated_padj_mean': self.downregulated_genes[padj_column].mean(),
                    'downregulated_padj_median': self.downregulated_genes[padj_column].median(),
                })
            
            self.summary_stats = stats
            logger.info("Summary statistics calculated successfully")
            return stats
        except Exception as e:
            logger.error(f"Error calculating summary statistics: {str(e)}")
            return None
    
    def display_summary_statistics(self) -> None:
        """Display summary statistics in a formatted manner."""
        if self.summary_stats is None:
            logger.warning("Summary statistics not calculated yet. Call calculate_summary_statistics() first.")
            return
        
        logger.info("=" * 80)
        logger.info("SUMMARY STATISTICS")
        logger.info("=" * 80)
        
        logger.info(f"Total genes analyzed: {self.summary_stats['total_genes']}")
        logger.info(f"Upregulated genes: {self.summary_stats['upregulated_count']} "
                   f"({self.summary_stats['upregulated_percent']:.2f}%)")
        logger.info(f"Downregulated genes: {self.summary_stats['downregulated_count']} "
                   f"({self.summary_stats['downregulated_percent']:.2f}%)")
        
        logger.info("\nUpregulated Genes - Log2FoldChange Statistics:")
        logger.info(f"  Mean: {self.summary_stats['upregulated_log2fc_mean']:.4f}")
        logger.info(f"  Median: {self.summary_stats['upregulated_log2fc_median']:.4f}")
        logger.info(f"  Std Dev: {self.summary_stats['upregulated_log2fc_std']:.4f}")
        logger.info(f"  Range: [{self.summary_stats['upregulated_log2fc_min']:.4f}, "
                   f"{self.summary_stats['upregulated_log2fc_max']:.4f}]")
        
        logger.info("\nDownregulated Genes - Log2FoldChange Statistics:")
        logger.info(f"  Mean: {self.summary_stats['downregulated_log2fc_mean']:.4f}")
        logger.info(f"  Median: {self.summary_stats['downregulated_log2fc_median']:.4f}")
        logger.info(f"  Std Dev: {self.summary_stats['downregulated_log2fc_std']:.4f}")
        logger.info(f"  Range: [{self.summary_stats['downregulated_log2fc_min']:.4f}, "
                   f"{self.summary_stats['downregulated_log2fc_max']:.4f}]")
        
        logger.info("=" * 80)
    
    def generate_volcano_plot(self, log2fc_column: str = 'Log2FoldChange',
                             padj_column: str = 'padj', log10_threshold: float = -np.log10(0.05),
                             fc_threshold: float = 1.0, figsize: Tuple[int, int] = (12, 8)) -> None:
        """
        Generate a volcano plot visualization.
        
        Args:
            log2fc_column (str): Name of Log2FoldChange column (default: 'Log2FoldChange')
            padj_column (str): Name of adjusted p-value column (default: 'padj')
            log10_threshold (float): -log10(p-value) threshold for significance (default: 1.3 for p=0.05)
            fc_threshold (float): Log2FoldChange threshold (default: 1.0)
            figsize (Tuple[int, int]): Figure size (default: (12, 8))
        """
        if self.data is None or padj_column not in self.data.columns:
            logger.error(f"Data or '{padj_column}' column not found")
            return
        
        try:
            plt.figure(figsize=figsize)
            
            # Calculate -log10(padj)
            data_plot = self.data.copy()
            data_plot['neg_log10_padj'] = -np.log10(data_plot[padj_column].replace(0, 1e-300))
            
            # Define colors
            colors = np.where(
                (data_plot[log2fc_column] > fc_threshold) & (data_plot['neg_log10_padj'] > log10_threshold),
                'red',  # Upregulated
                np.where(
                    (data_plot[log2fc_column] < -fc_threshold) & (data_plot['neg_log10_padj'] > log10_threshold),
                    'blue',  # Downregulated
                    'gray'  # Not significant
                )
            )
            
            plt.scatter(data_plot[log2fc_column], data_plot['neg_log10_padj'], 
                       c=colors, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
            
            # Add threshold lines
            plt.axhline(y=log10_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
            plt.axvline(x=fc_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
            plt.axvline(x=-fc_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
            
            plt.xlabel(f'{log2fc_column}', fontsize=12, fontweight='bold')
            plt.ylabel(f'-log10({padj_column})', fontsize=12, fontweight='bold')
            plt.title('Volcano Plot: Gene Expression Changes', fontsize=14, fontweight='bold')
            plt.grid(True, alpha=0.3)
            
            # Add legend
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='red', label='Upregulated'),
                Patch(facecolor='blue', label='Downregulated'),
                Patch(facecolor='gray', label='Not Significant')
            ]
            plt.legend(handles=legend_elements, loc='upper right')
            
            output_path = self.output_dir / 'volcano_plot.png'
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Volcano plot saved to {output_path}")
            plt.close()
        except Exception as e:
            logger.error(f"Error generating volcano plot: {str(e)}")
    
    def generate_bar_charts(self, n: int = 15, log2fc_column: str = 'Log2FoldChange',
                           figsize: Tuple[int, int] = (16, 10)) -> None:
        """
        Generate bar charts for top upregulated and downregulated genes.
        
        Args:
            n (int): Number of top genes to display (default: 15)
            log2fc_column (str): Name of Log2FoldChange column (default: 'Log2FoldChange')
            figsize (Tuple[int, int]): Figure size (default: (16, 10))
        """
        if self.upregulated_genes is None or self.downregulated_genes is None:
            logger.error("Genes not separated yet")
            return
        
        try:
            fig, axes = plt.subplots(2, 1, figsize=figsize)
            
            # Get gene names
            def get_gene_name(row):
                return row.get('Gene', row.get('gene_id', str(row.name)))
            
            # Top upregulated genes
            top_up = self.upregulated_genes.nlargest(n, log2fc_column).copy()
            top_up['gene_name'] = top_up.apply(get_gene_name, axis=1)
            top_up = top_up.sort_values(log2fc_column)
            
            axes[0].barh(range(len(top_up)), top_up[log2fc_column], color='darkred', alpha=0.7)
            axes[0].set_yticks(range(len(top_up)))
            axes[0].set_yticklabels(top_up['gene_name'], fontsize=9)
            axes[0].set_xlabel(f'{log2fc_column}', fontsize=11, fontweight='bold')
            axes[0].set_title(f'Top {n} Upregulated Genes', fontsize=12, fontweight='bold')
            axes[0].grid(axis='x', alpha=0.3)
            
            # Top downregulated genes
            top_down = self.downregulated_genes.nsmallest(n, log2fc_column).copy()
            top_down['gene_name'] = top_down.apply(get_gene_name, axis=1)
            top_down = top_down.sort_values(log2fc_column, ascending=False)
            
            axes[1].barh(range(len(top_down)), top_down[log2fc_column], color='darkblue', alpha=0.7)
            axes[1].set_yticks(range(len(top_down)))
            axes[1].set_yticklabels(top_down['gene_name'], fontsize=9)
            axes[1].set_xlabel(f'{log2fc_column}', fontsize=11, fontweight='bold')
            axes[1].set_title(f'Top {n} Downregulated Genes', fontsize=12, fontweight='bold')
            axes[1].grid(axis='x', alpha=0.3)
            
            plt.tight_layout()
            output_path = self.output_dir / 'top_genes_bar_charts.png'
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Bar charts saved to {output_path}")
            plt.close()
        except Exception as e:
            logger.error(f"Error generating bar charts: {str(e)}")
    
    def generate_distribution_plot(self, log2fc_column: str = 'Log2FoldChange',
                                  figsize: Tuple[int, int] = (14, 6)) -> None:
        """
        Generate distribution plots for Log2FoldChange values.
        
        Args:
            log2fc_column (str): Name of Log2FoldChange column (default: 'Log2FoldChange')
            figsize (Tuple[int, int]): Figure size (default: (14, 6))
        """
        if self.data is None:
            logger.error("Data not loaded")
            return
        
        try:
            fig, axes = plt.subplots(1, 2, figsize=figsize)
            
            # Histogram
            axes[0].hist(self.data[log2fc_column], bins=50, color='steelblue', alpha=0.7, edgecolor='black')
            axes[0].axvline(x=0, color='red', linestyle='--', linewidth=2, label='No change')
            axes[0].set_xlabel(f'{log2fc_column}', fontsize=11, fontweight='bold')
            axes[0].set_ylabel('Frequency', fontsize=11, fontweight='bold')
            axes[0].set_title(f'Distribution of {log2fc_column}', fontsize=12, fontweight='bold')
            axes[0].legend()
            axes[0].grid(alpha=0.3)
            
            # Density plot
            self.data[log2fc_column].plot(kind='density', ax=axes[1], color='steelblue', linewidth=2)
            axes[1].axvline(x=0, color='red', linestyle='--', linewidth=2, label='No change')
            axes[1].fill_between(axes[1].get_lines()[0].get_xdata(), 
                                axes[1].get_lines()[0].get_ydata(), alpha=0.3, color='steelblue')
            axes[1].set_xlabel(f'{log2fc_column}', fontsize=11, fontweight='bold')
            axes[1].set_ylabel('Density', fontsize=11, fontweight='bold')
            axes[1].set_title(f'Density Plot of {log2fc_column}', fontsize=12, fontweight='bold')
            axes[1].legend()
            axes[1].grid(alpha=0.3)
            
            plt.tight_layout()
            output_path = self.output_dir / 'distribution_plots.png'
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Distribution plots saved to {output_path}")
            plt.close()
        except Exception as e:
            logger.error(f"Error generating distribution plots: {str(e)}")
    
    def export_results(self, log2fc_column: str = 'Log2FoldChange') -> None:
        """
        Export analysis results to CSV files.
        
        Args:
            log2fc_column (str): Name of Log2FoldChange column (default: 'Log2FoldChange')
        """
        try:
            # Export upregulated genes
            if self.upregulated_genes is not None:
                up_path = self.output_dir / 'upregulated_genes.csv'
                upregulated_sorted = self.upregulated_genes.sort_values(log2fc_column, ascending=False)
                upregulated_sorted.to_csv(up_path, index=False)
                logger.info(f"Upregulated genes exported to {up_path}")
            
            # Export downregulated genes
            if self.downregulated_genes is not None:
                down_path = self.output_dir / 'downregulated_genes.csv'
                downregulated_sorted = self.downregulated_genes.sort_values(log2fc_column, ascending=True)
                downregulated_sorted.to_csv(down_path, index=False)
                logger.info(f"Downregulated genes exported to {down_path}")
            
            # Export full results with classification
            if self.data is not None:
                full_path = self.output_dir / 'all_genes_classified.csv'
                classified_data = self.data.copy()
                classified_data['regulation'] = classified_data[log2fc_column].apply(
                    lambda x: 'upregulated' if x > 0 else ('downregulated' if x < 0 else 'unchanged')
                )
                classified_data = classified_data.sort_values(log2fc_column, ascending=False)
                classified_data.to_csv(full_path, index=False)
                logger.info(f"Classified genes exported to {full_path}")
            
            # Export summary statistics
            if self.summary_stats is not None:
                stats_path = self.output_dir / 'summary_statistics.csv'
                stats_df = pd.DataFrame(list(self.summary_stats.items()), columns=['Statistic', 'Value'])
                stats_df.to_csv(stats_path, index=False)
                logger.info(f"Summary statistics exported to {stats_path}")
        except Exception as e:
            logger.error(f"Error exporting results: {str(e)}")
    
    def run_full_analysis(self, input_file: str = None, log2fc_column: str = 'Log2FoldChange',
                         padj_column: str = 'padj') -> bool:
        """
        Run the complete gene regulation analysis pipeline.
        
        Args:
            input_file (str): Path to input Excel file (optional, uses initialized path if not provided)
            log2fc_column (str): Name of Log2FoldChange column (default: 'Log2FoldChange')
            padj_column (str): Name of adjusted p-value column (default: 'padj')
        
        Returns:
            bool: True if successful, False otherwise
        """
        if input_file:
            self.input_file = input_file
        
        logger.info("Starting full gene regulation analysis...")
        logger.info("=" * 80)
        
        # Step 1: Load data
        if not self.load_data():
            return False
        
        # Step 2: Separate genes
        if not self.separate_genes(log2fc_column):
            return False
        
        # Step 3: Calculate statistics
        self.calculate_summary_statistics(log2fc_column, padj_column)
        self.display_summary_statistics()
        
        # Step 4: Display top genes
        self.display_top_genes(n=10, log2fc_column=log2fc_column, padj_column=padj_column)
        
        # Step 5: Generate visualizations
        logger.info("Generating visualizations...")
        self.generate_volcano_plot(log2fc_column=log2fc_column, padj_column=padj_column)
        self.generate_bar_charts(n=15, log2fc_column=log2fc_column)
        self.generate_distribution_plot(log2fc_column=log2fc_column)
        
        # Step 6: Export results
        logger.info("Exporting results...")
        self.export_results(log2fc_column=log2fc_column)
        
        logger.info("=" * 80)
        logger.info("Gene regulation analysis completed successfully!")
        logger.info(f"All results saved to: {self.output_dir.absolute()}")
        
        return True


def main():
    """Main execution function."""
    # Define input and output paths
    input_file = 'TFR-OE--Control_sig-genes.xlsx'
    output_directory = 'results'
    
    # Create analyzer instance
    analyzer = GeneRegulationAnalyzer(input_file=input_file, output_dir=output_directory)
    
    # Run full analysis
    success = analyzer.run_full_analysis(
        log2fc_column='Log2FoldChange',
        padj_column='padj'
    )
    
    if not success:
        logger.error("Analysis failed. Please check the input file and logs.")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
