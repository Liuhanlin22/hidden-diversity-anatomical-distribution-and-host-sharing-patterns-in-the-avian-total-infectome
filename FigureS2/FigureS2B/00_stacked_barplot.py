import pandas as pd
import matplotlib
# Use a non-interactive backend for batch rendering.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# ========================================================
# Keep text editable in vector output where supported.
# ========================================================
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
# ========================================================

def main():
    input_file = "00_qc.xlsx"
    output_file = "Family_Composition_by_Order.pdf"
    
    print(f"Reading input data: {input_file} ...")
    
    try:
        df = pd.read_excel(input_file)
    except FileNotFoundError:
        print(f"Error: input file not found: '{input_file}'.")
        return

    # Remove accidental whitespace from column names.
    df.columns = df.columns.str.strip()

    # Expected column names in the input workbook.
    col_order = 'Order'
    col_family = 'Family'

    if col_order not in df.columns or col_family not in df.columns:
        print(f"Error: missing '{col_order}' or '{col_family}' in the input table.")
        print(f"Available columns: {df.columns.tolist()}")
        sys.exit(1)

    print("Calculating Family composition within each Order ...")
    
    # Count rows within each Order and normalize to percentages.
    # To weight by an abundance column, replace this with a weighted crosstab.
    
    cross_tab = pd.crosstab(df[col_order], df[col_family], normalize='index') * 100

    print("Rendering the 100% stacked bar plot ...")
    
    # 1. Match the dimensions used by the R implementation.
    fig, ax = plt.subplots(figsize=(20, 7))

    # 2. Generate a distinct palette for the Family categories.
    num_families = len(cross_tab.columns)
    colors = sns.color_palette("husl", num_families).as_hex()

    # 3. Draw the stacked bars.
    cross_tab.plot(kind='bar', stacked=True, ax=ax, color=colors, width=0.7, edgecolor='none')

    # 4. Apply the plot theme.
    ax.set_facecolor('white') # panel.background = element_rect(fill = "white")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('#333333')
    ax.spines['left'].set_color('#333333')
    
    ax.set_title('', fontweight='bold')
    ax.set_xlabel('Order', fontweight='bold', fontsize=12)
    ax.set_ylabel('Percentage (%)', fontweight='bold', fontsize=12)

    # axis.text.x = element_text(angle = 45, hjust = 1)
    plt.xticks(rotation=45, ha='right')

    # legend.position = "right", legend.text = element_text(size = 8)
    # Place the legend outside the plotting area.
    ax.legend(title='Family', bbox_to_anchor=(1.01, 1), loc='upper left', 
              fontsize=8, frameon=False, 
              # Use two columns for a large number of Family categories.
              ncol=1 if num_families <= 25 else 2) 

    # Keep long legends and axis labels within the exported canvas.
    plt.tight_layout()

    # 5. Save the figure.
    plt.savefig(output_file, dpi=300, format='pdf', bbox_inches='tight', transparent=True)
    
    print("Plot completed.")
    print(f"Figure saved to: {output_file}")

if __name__ == "__main__":
    main()
