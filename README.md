# Project Proposal: **Mapping and Understanding the Fight Against Cancer Using BindingDB**
## Abstract (150 words)
This project explores cancer treatment evolution through an analysis of cancer-related data from BindingDB, focusing on how research trends have changed over time and what these trends reveal about treatment progress. This project aims at understanding relationships between ligands and proteins to extract valuable information on cancer treatment. By investigating how the proportions of treatments vary among cancer types, as well as the impact of treatments targeting mutant proteins, we assess how current approaches address challenges posed by cancer mutations, pathways, and side effects. Our goal is to tell the story of cancer treatment advancements, identifying pivotal research milestones and revealing insights on the struggles due to side effects. This work ultimately illustrates the role of molecular interaction data in advancing cancer research and addressing critical health challenges.

## Research Questions 	 
1. How has the number of treatments evolved through time? Miki and Kam
2. Are there any cancer-related proteins that rapidly became major research targets? Miki and Kam 
3. Which cancer-related proteins are most frequently targeted by therapeutic drugs, based on BindingDB and DrugBank data? Mathis
4. How many mutants does every cancer-related protein have, and how effective are the ligands against these mutants? Matheo 
5. What are the key ligands that ultimately get approved as drugs, and what are there binding/inhibiting/concatring/kinetic properties? Seb 
6. What are the side effects of ligands that show great properties against cancer ? Seb

## Additional datasets
While **BindingDB** will serve as the primary dataset, **DrugBank** provides valuable information on treatment approval, side effects, and other pharmacological properties of the ligands. 
### 1. **DrugBank** (Additional Dataset)
- **Source**: DrugBank (https://www.drugbank.ca/)
- **Size and Format**: Structured data (XML) with detailed drug information, including targets, mechanisms of action, and clinical uses.
- **Description**: DrugBank contains extensive data on over 13,000 drugs and their interactions with biological targets. Integrating **DrugBank** data with **BindingDB** allows us to explore how specific drugs interact with cancer-related proteins, aiding in the identification of potential therapeutic compounds. Additionally, Drugbank contains the information on approved drugs and hence, ligands that were deemed potent enough, as well their side effects.
- **Plan for Management**: Data from **DrugBank** will be mapped to relevant proteins from **BindingDB** based on their interaction data (unique names, or structure), allowing us to identify high-affinity binding drugs and evaluate their effectiveness in treating cancer.
**Additional sources for context and story-telling** : SEER, COSMIC, TCGA, IARC

## Methods
### 1. **Data Preprocessing**:
- Clean and preprocess data from **BindingDB**, **DrugBank**,
-  Map Protein IDs in **BindingDB** to those of **DrugBank** to integrate drug-target interactions contained in BindingDB to clinical data from DrugBank. IDs used : 
1. PubChem CID: PubChem, CID, Compound, 
2. ChEBI ID of Ligand: ChEBI, Ligand, Bioentity,
3. ChEMBL ID of Ligand: ChEMBL, Ligand, Bioactivity,
4. DrugBank ID of Ligand: DrugBank, Drug, Pharmaceutical,
5. KEGG ID of Ligand: KEGG, Ligand, Pathway,
6. ZINC ID of Ligand: ZINC, Ligand, Screening, 
7. Ligand SMILES: SMILES, Structure, Notation
8. Ligand InChI Key: InChI Key, Structure, Identifier
9. BindingDB MonomerID: BindingDB, Monomer, Interaction,
- Merge **BindingDB**, **DrugBank** to and save the newly formed dataset : Merged_DB . 
- Filter the Merged_DB based on its “specific-function”, “Target Name”, “entry name of target chain” attributes (for proteins) by identifying cancer-related keywords. 

### 2. **Exploratory Data Analysis (EDA)**:
- Look into correlations between all attributes of interest, in particular binding properties, temperature and pH. 
- Visualize the distribution of binding affinities for cancer-related proteins and the drugs associated with them using histograms.

### 3 **Protein Analysis**:
- Perform trend analysis by looking at the proteins present in published articles using their DOI.

### 4. **Network Analysis**:
- Link proteins to their mutants based on their target names and find all ligands that were tested on these proteins and have properties above a certain threshold. 

### 5. **Ligand Analysis**:
- Compare chemical properties of ligands, and identify which ones are present in approved medical drugs. 
- Analyze their side effects, by looking at the other molecules they bind to. 

### 6. **Visualization**:
- Create interactive visualizations using tools like **Plotly**, **Seaborn**, or **Matplotlib** to display protein-ligand interactions, treatment efficacy, and survival outcomes.
- Use network diagrams to represent **protein-drug interactions** (4).
## Proposed timeline
Week 10 : Identify which hypotheses are true, while analyzing time trends, and proteins (visualization)
Week 11 :  Network Analysis, Ligand Analysis (visualization)
Week 12 :  Side effects (visualization)
Week 13 : Conclude and present the storyline on our website.
## Organization within the team
Each of us is responsible for one topic. Each week, we talk about our findings and present our visualizations.
## Completed milestones	 	
Find data on cancer cases, deadliness over the years, and trends in general
Get access to DrugBank and understand how to link it to BindingDB 
Create template to read BindingDB efficiently > mySQL > tsv 
Preprocessing of the data - data cleaning 	
Identify proteins, which are related to cancer
Get a better understanding on Chemical cancer treatment to understand the main metabolic pathways 
Define all cancer-related parameters that could show interesting trends over the years and will help us lay out the story of the fight against cancer 

