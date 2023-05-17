import streamlit as st
import pandas as pd
import plotly.express as px

def load_data():
    df = pd.read_csv('raw_data.txt', delimiter='\t', encoding='ISO-8859-1', skiprows=27)
    df = df[df['Name'].notna()] 
    df = df[~df['Name'].isin(['Blank', 'Empty'])] 
    df = df[df['ID'].notna()] 
    df = df[df['ID'] != 'EMPTY']  
    df['F635 Corrected'] = df['F635 Mean'] - df['B635 Mean']
    df['F532 Corrected'] = df['F532 Mean'] - df['B532 Mean']
    df['F488 Corrected'] = df['F488 Mean'] - df['B488 Mean']
    df['Signal'] = df[['F635 Corrected', 'F532 Corrected', 'F488 Corrected']].mean(axis=1)
    control_df = df[df['Name'] == 'Control_1']
    control_signal_mean = control_df['Signal'].mean()
    df['Signal'] = df['Signal'] - control_signal_mean
    all_df = df.copy()  # Store a copy of the data before filtering
    df = df[df['Signal'] > 0]
    return df, all_df

df, all_df = load_data()

st.title("Epitope Mapping App")

quantile_value = st.slider('Select Threshold Percentile:', min_value=0.30, max_value=1.0, value=0.80, step=0.05)

df_mean = df.groupby('ID')['Signal'].mean().reset_index()
threshold = df_mean['Signal'].quantile(quantile_value)
df_mean['Positive'] = df_mean['Signal'] > threshold
positive_peptides_ids = df_mean[df_mean['Positive']]['ID'].tolist()
positive_df = df[df['ID'].isin(positive_peptides_ids)]

# First plot with only positive peptides, color gradient is red with higher signal
st.markdown("## Plot of Positive Peptides")
fig1 = px.scatter(positive_df, x='Column', y='Row', color='Signal', color_continuous_scale='Reds', hover_data={'Signal': ':.0f', 'Name': True, 'ID': True})
st.plotly_chart(fig1)



threshold_F635 = df['F635 Corrected'].quantile(quantile_value)
threshold_F532 = df['F532 Corrected'].quantile(quantile_value)
threshold_F488 = df['F488 Corrected'].quantile(quantile_value)

nanobody_epitopes = df.loc[df['F488 Corrected'] > threshold_F488, 'ID'].unique()
serum_epitopes = df.loc[df['F635 Corrected'] > threshold_F635, 'ID'].unique()
common_epitopes = set(nanobody_epitopes).intersection(serum_epitopes)

common_epitopes_df = df[df['ID'].isin(common_epitopes)]
common_epitopes_df = common_epitopes_df.groupby('ID').agg({'F635 Corrected': 'mean', 'F488 Corrected': 'mean'}).reset_index()
common_epitopes_df = common_epitopes_df.round(0) # Remove decimals
common_epitopes_df['F635 Corrected'] = common_epitopes_df['F635 Corrected'].astype(int)
common_epitopes_df['F488 Corrected'] = common_epitopes_df['F488 Corrected'].astype(int)
common_epitopes_df = common_epitopes_df.sort_values(by=['F635 Corrected', 'F488 Corrected'], ascending=False)
# common_epitopes_df[['F635 Corrected', 'F488 Corrected']] = common_epitopes_df[['F635 Corrected', 'F488 Corrected']].format("{:.0f}")



st.markdown("## Table of Common Epitopes between nanobody and serum")



#Download button for the table in FASTA file format
fasta_data = ""
for seq_id in common_epitopes_df['ID']:
    fasta_data += f">{seq_id}\n"
    fasta_data += f"{seq_id}\n"


st.download_button(
label="Download FASTA",
data=fasta_data,
file_name='common_epitopes.fasta',
mime='text/plain'
)

# CSS to inject contained in a string
hide_table_row_index = """
            <style>
            thead tr th:first-child {display:none}
            tbody th {display:none}
            </style>
            """
# Inject CSS with Markdown
st.markdown(hide_table_row_index, unsafe_allow_html=True)

st.table(common_epitopes_df)

# Second plot with all data, color gradient is red with higher signal
st.markdown("## Plot of All Data")
fig2 = px.scatter(all_df, x='Column', y='Row', color='Name', hover_data=['Name', 'ID'] )
st.plotly_chart(fig2)