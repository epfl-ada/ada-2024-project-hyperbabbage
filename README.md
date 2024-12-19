# Project Proposal: **Mapping and Understanding the Fight Against Cancer Using BindingDB**

Please read our [data story](http://217.160.247.121/).

## Abstract
This project explores cancer treatment evolution through an analysis of cancer-related data from BindingDB, focusing on how research trends have changed over time and what these trends reveal about treatment progress. This project aims at understanding relationships between ligands and proteins to extract valuable information on cancer treatment. By investigating how the proportions of treatments vary among cancer types, as well as the impact of treatments targeting mutant proteins, we assess how current approaches address challenges posed by cancer mutations, pathways, and side effects. Our goal is to tell the story of cancer treatment advancements, identifying pivotal research milestones and revealing insights on the struggles due to side effects. This work ultimately illustrates the role of molecular interaction data in advancing cancer research and addressing critical health challenges.

## Research Questions 	 
1. How has the number of treatments evolved through time?
2.  Are there any cancer-related proteins that rapidly became major research targets?
3.  Which cancer-related proteins are most frequently targeted by therapeutic drugs, based on BindingDB and DrugBank data?
4. How many mutants does every cancer-related protein have, and how effective are the ligands against these mutants? 
5. What are the key ligands that ultimately get approved as drugs, and what are there binding/inhibiting/concatring/kinetic properties? 
6. How do mutations in cancer-related proteins influence the binding affinity of therapeutic drugs and treatment success, knowing that cancer cells mutate to gain treatment resistance?
7. What are the side effects of ligands that show great properties against cancer ? 

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
## Contributions
  - Mikuláš Vanoušek:
    - Preprocessing BindingDB
    - Collecting metadata about research
    - Time trends analysis
    - GPT for cancer-related proteins identification
    - Content of the data story
  - Sebastian Delsad:
    - Merging Binding DB with DrugBank
    - Preprocessing the merged dataset
    - Network analysis
    - Code for Cancer Related Proteins
    - The data story website
  - Kamel Charaf
    - Preprocessing BindingDB
    - Collecting metadata about research
    - Time trends analysis
    - GPT for cancer-related proteins identification
    - Massive code cleanup
  - Mathis Finckh : 
    - Initial trials on cancer-related protein search : 
    - through related cancer-drugs (web-scrapping)
    - through keyword search
    - Toxicity analysis
    - Part of the content of the data story
  - Matheo Godenzi : 
    * storyline idea 
    * ⁠collecting metadata 
    * ⁠merge review
    * ⁠chemical analysis 
    * ⁠cancer process identification and analysis 
    * ⁠website content writing



## Running the analysis
All results are produced by the `./results.ipynb` notebook.

### Enviornment
Create a the conda environment by running the following command in the root of the project:
```bash
conda env create --prefix ./.conda -f ./environment.yml
```

Activate the environment by running:
```bash
conda activate ./.conda
```

### Data
You can download both clean and processed data from [Google Drive](https://drive.google.com/file/d/1J1f1xcV4c9FIRCuid7s0Gyv6gxitS68N/view?usp=drive_link). Please place the `data` folder in the root directory of the project.

If you would like to create the clean data yourself, delete the contents of `data/clean` folder and run `./preprocessing.ipynb`. Please note that you will need ~32GB of memory to merge BindingDB and DrugBank, and fetching the metadata about research based on DOI from an API takes at least an hour.
