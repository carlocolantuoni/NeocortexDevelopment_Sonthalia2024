# NeMO-Analytics
This document introduces the NeMO (Neuroscience Multi-Omic Archive) environment workflow and highlights the main features of the NeMO Analytics platform, providing a quick overview for new users. NeMO Analytics is closely related to the Gene Expression Analysis Resource (gEAR), as both platforms are built on the same back-end data infrastructure. The key distinction lies under the scope of the study, where NeMO is designed specifically for neuroscience study, and gEAR is primarily tailored to hearing-related study. Because they share the same infrastructure, many of the features and workflows are common across both platforms.

> *This tutorial is adapted from the official gEAR wiki page (which can be accessed through this [link](https://github.com/IGS/gEAR/blob/main/docs/Documentation/gEARWiki.md)), with modifications to highlight the NeMO-specific environment and workflows.*

---

## Overview of typical NeMO workflow 
<img width="715" height="627" alt="image" src="https://github.com/user-attachments/assets/422d7fdc-26b5-4536-98c6-d461316302fb" />

---

### 1. Getting started with NeMO platform 

#### 1.1 Set up NeMO accounts
- Go to the homepage → click the drop-down button next to "Log In" on the top right → select Sign Up.
- Enter your basic details and email address.
- Verify your email to activate the account.
- Once logged in, you will have access to all the NeMO analytic tools and be able to make private/public dataset collections.
- Manage your profile anytime by clicking your username (top right).
- Enable colorblind-friendly mode for plots if needed.

<img width="1378" height="654" alt="image" src="https://github.com/user-attachments/assets/0deee0e6-fdc4-4c0a-a5f3-43b01d4d2fb7" />

---

#### 1.2 Layout of the NeMO platform 
    
*1.2.1 Home page* 

The homepage of the NeMO platform is the "My Workspace" Dashboard. At any time, you can return to the homepage by clicking the NeMO icon in the top left corner. 

<img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/1695c719-5d58-41c7-8a4c-ca0a28a5e7a9" />

*1.2.2 Navigation panel*

The sidebar on the left-hand side of the page is the navigation panel. This panel will stay on the side as the users move between pages, allowing quick access to the various modules within NeMO. As the users become more familiar with the layout of the panel, they can choose to minimize it so that only the symbols are displayed (without text labels). To do this, simply click the arrow located in the upper-right corner of the navigation panel.

---

#### 1.3 Explore and exploit the available datasets

As a comprehensive data compendium, the NeMO platform allows users not only to upload their own datasets for customized analysis but also to explore a wide range of publicly available datasets. In fact, one of the best ways to become familiar with NeMO environment is to dive into the extensive collection of datasets already deposited on the platform. This section aims to introduce the core functions and steps for data exploration, helping new users build an understanding of the platform and practice analyzing datasets of interest without uploading their own data.

*1.3.1 Finding datasets*

To view results from a specific dataset or browse available datasets, you can click the “Datasets” option or "Dataset collections" from the navigation panel. Within the Dataset Explorer, you can refine searches using the filter controls on the left side of the page. Filters include options such as data type, keywords from the study title, or organism. The explorer supports four different views of datasets, including "Table View", "List View-Compact", "List View (Expanded), and "Collection View". Each view presents a different level of detail about the datasets. 

<img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/ab54a4eb-9831-4baf-a83f-e859e65c01f5" />

To search for specific genes within a dataset,
  - Click the small arrow in front of the dataset title to expand the entry.
  - In the expanded entry, you will see a default display created by the dataset owner.
  - Click the test tube icon to open the Expression Viewer, which allows you to explore gene expression using the default display.

If you prefer to create your own visualization:
  - Use the single-gene curator option from the drop-down menu under the analysis tools.
  - Detailed instructions on building custom visualizations will be provided in "5. Curate Visualization"

When you first open the Expression Viewer, a default gene will appear, but you can replace it by typing any other gene name into the search box. By default, the viewer shows one gene at a time with the "single-gene display" mode. However, after saving your own gene list in NeMO (covered in a later section), you will be able to view the saved list of genes  using the "multi-gene display" mode if multi-gene displays have been saved in the selected collection or dataset. 

<img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/357a677a-d075-4052-8b01-2a927b33a5ac" />

<img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/c2931cc9-002c-460a-ae84-c1c193a6fc1c" />

On the right side of the Expression Viewer display, there are two buttons:
  - The left button lets you zoom in on the visualization.
  - The right button opens additional tools, including the option to view or switch to other displays created by the dataset owner.
    - Click on the "Choose display" button allows you to select the new display you wish to view your gene of interest.
    - Other buttons are also available for users to interact with displays and datasets using other NeMO analytic tools.  
      
<img width="1378" height="452" alt="image" src="https://github.com/user-attachments/assets/f6da7f8a-1591-4c29-af8a-ca43204261af" />


---

### 2. Retrieve datasets and metadata
Step 2 of the NeMO workflow is mostly conducted offline, where users can identify their own interested datasets and acquire processed gene-level expression data and exhaustive metadata information. The next few sections will begin introducing how to upload your own multi-omics datasets. Step 2 of the NeMO workflow is mostly user-dependent and conducted offline. 

---

### 3. Format standardizing and NeMO output packaging 
`Counts.IN.SeuratDimRedNeMO.OUT.TRIM.R` is a R function we have implemented in the NeMO workflow to standardize the multi-omics datasets format and package them into `NeMOanalytics.org` compatible outputs for uploading. Details about how to use the function can be found in `Counts.IN.SeuratDimRedNeMO.OUT.TRIM.md`. 

---

### 4. Upload output to NeMO analytics platform

**4.1 Output from the NeMO function**

Running the NeMO function will generate four outputs in the designated directory. Download the .tar.gz file and the .xlsx for uploading it onto NEMO

<img width="450" height="122" alt="image" src="https://github.com/user-attachments/assets/5380585e-f975-4cc4-a0a6-92bc3b4ca119" />

---

**4.2 Uploading Data to NeMO**

The dataset upload process in NeMO consists of two main parts: uploading the metadata file and uploading the expression data file.

*4.2.1. Upload Metadata*
 - In the left navigation panel, click Dataset Uploader.
 - On the Dataset Upload page, click “Choose” and select the metadata .xlsx file generated by the NeMO function.
 - Confirm your selection by clicking Upload.
 - The metadata fields will automatically populate, except for the Dataset type.
 - Manually select the correct dataset type from the dropdown menu.
 - Click Submit to proceed.
   <img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/632d158e-7670-45d0-9258-f747e97c7d4a" />
---

*4.2.2 Upload Expression Data*

- On the next page, under Step – Upload dataset, select MEX / 3-tab format.
- Click “Choose” and upload the .tar.gz file generated by the NeMO function.
- Click Upload dataset to begin.
- The upload time will depend on file size. Do not close your browser until the upload finishes, as this will interrupt the process.
<img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/b97826cb-c4a3-4047-853e-f569721ab9e0" />

*4.2.3. Background Processing*
   - Once the dataset file is successfully uploaded, NeMO will automatically prompt you to the next step, "processing", where all the extraction and processing of the dataset happens in the background. At this stage, you no longer need to keep your browser open.
   - When processing is complete, click continue to move on to the next setp.
<img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/a2e290c9-5c0e-41f6-a0db-44ce97df46ec" />

*4.2.4 Finalize and Curate*

- After processing is complete, follow the prompt on NeMO to change the data visibility and finalize the uploading. 
- You can then begin curating dataset visualizations and adding the dataset to your collections.

<img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/25ba3da4-6804-4369-9788-bd0588cfc242" />
<img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/287b9247-43df-4d3a-ba19-f7810e6dc525" />


*4.2.5 Troubleshooting*

If you encounter errors during the upload process, please submit a ticket on the gEAR GitHub site
(https://github.com/IGS/gEAR/issues). Be sure to include your dataset ID when reporting the issue.

---

### 5. Curate Visualization 

One of the major strengths of the NeMO Analytics platform is the ability to create multiple visualizations for the same dataset. In addition to the default displays generated by the original dataset owner (as shown in 1.3 Explore and exploit the available datasets), you can create custom displays that are better suited to your own research questions for the datasets that are available on NeMO or your own uploaded datasets. 

**5.1 Launch the visualization tools**

There are three ways you could start curating your own visualization, depending on the usage scenario
 1. After uploading your own datasets

Clicking on the "Curate" button will bring you to the "single-gene display" page, where your uploaded dataset will automatically populate in the "1. select a dataset" slot

<img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/674640b7-b337-40cd-aa49-4b0d91090f99" />

 2. Customize visualization for existing datasets

-  The "single-gene display" page can also be accessed through the navigation panel on the left. When the "single-gene display" is launched, users can  manually search for the dataset of interest from the drop-down menu or type in the name of the dataset.
<img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/25ea897d-f1a9-4735-97dd-f4a8f6bf2ed5" />

- Alternatively, in the "Dataset collection" page, select the "single-gene curator" in the expanded entry of any dataset from the "analysis tool" dropdown menu. This will also bring the user to the "single-gene display" page with the name of the selected dataset populated. 

<img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/1fdd31af-844f-4cfb-a69c-c8fe10cfceca" />

---

**5.2 Create new visualization**
After navigating to the "single-gene display" page, you can start creating a new visualization for the dataset by clicking "Curate new display".  

<img width="1378" height="443" alt="image" src="https://github.com/user-attachments/assets/8dc5542a-4b6a-408f-a91f-57257bc50e6e" />


Multiple plot types and analysis types are available through NeMO, including line plot, bar plot, scatter plot, violin plot, t-SNE, UMAP, and PCA. NeMO also offers the option to make the plot interactive or static. All these features for the plot can be changed in step three, "Select dataset plot type and analysis type". 


<img width="1378" height="184" alt="image" src="https://github.com/user-attachments/assets/355b4d98-8f5a-4ec9-be14-a07abee9fd8d" />

Choose a gene for your analysis in Step 4. Note that this gene could be arbitrary for the purpose of creating displays. Saving the visualization would allow you to look up a different gene/list of genes in the "Gene Expression Search" page (see 1.3 Explore and exploit the available datasets for details). 

The plot setting could be further customized by specifying the axes of the plots. The axis and color options available in the dropdown menu are related to the column names/variables in the cellMETA when running the `Counts.IN.SeuratDimRedNeMO.OUT.TRIM.R function`. Therefore, we recommend concatenating the complete metadata information to the cellMETA when first running the `Counts.IN.SeuratDimRedNeMO.OUT.TRIM.R function` to avoid going back and forth with the uploading process.  

<img width="1378" height="298" alt="image" src="https://github.com/user-attachments/assets/13d86a28-bd6f-411e-bfcf-5ee95a8ec9ea" />

Confirm the setting by clicking "Plot", which will bring you to the preview of the display, where you have more options to refine and make minor adjustments on the plot using "Dataset filters", "Plot configuration", and "Change color". Check the "Make this my default display box" if you wish it to appear as a preview on the expanded dataset list view. Confirm the change by clicking the "Save as new display" button. 
<img width="1378" height="525" alt="image" src="https://github.com/user-attachments/assets/28d50b78-73e8-4116-822e-398a1886905f" />


Once the display is saved, it will appear in the "single-gene display" page. To make changes to this existing display, select "Clone" under the display. After the necessary changes on the plot, confirm by selecting the "Update plot" button and choosing "Overwrite existing display instead". 
<img width="1378" height="526" alt="image" src="https://github.com/user-attachments/assets/552473ec-125b-420e-b2ec-72cc35708c00" />

---

### 6. Create dataset collections
Dataset collection is one of the unique NeMO features of grouping different datasets around a shared theme. The NeMO community has been curating dataset collections that are available to the public domain, but any user can create their own dataset collection of datasets to explore.

**6.1 Exploring site curated collections**

In the Gene Expression Viewer, you can search within specific dataset collections by using the dropdown menu. The menu provides several options, including: site-curated collection, collections you own, collections you’ve recently accessed, collections from groups you belong to, and collections shared with you by other users. 

Type in your gene of interest and simply select the desired collection from the dropdown to focus your search within that set of datasets.

<img width="1378" height="648" alt="image" src="https://github.com/user-attachments/assets/98c5f9b5-a53d-4a92-8ef0-ad6574f29c75" />

To generate a shareable link for this particular collection page, select Datasets or "Dataset collections" from the navigation panel. Use the dropdown menu from the Collection management to choose the desired curated collection. Move the cursor on the middle icon to get the shareable link. 
<img width="1374" height="648" alt="image" src="https://github.com/user-attachments/assets/78922e20-1999-4899-8afd-a3669dc0d869" />

---
**6.2  Create a new collection**

In the "Dataset Explorer" page, locate the "Collection Management" and click the plus (+) icon to create a new dataset collection. Enter a name for your collection and click "Add". The new collection will then appear in the dropdown menu.

<img width="1374" height="648" alt="image" src="https://github.com/user-attachments/assets/d5baed23-a9cb-4d56-96f0-0b8de6615fd7" />

---

**6.3 Add datasets to your collection**

In the same "Dataset Explorer" page, make sure you are in the collection you just created. Find the dataset you want to add to the collection and expand its dataset entry. Clicking the plus (+) icon next to a dataset entry adds it to the collection. To remove a dataset, click the minus (–) icon. 
<img width="1374" height="648" alt="image" src="https://github.com/user-attachments/assets/86fa18c3-a950-4306-b9a3-2c230bb064b8" />

Clicking the plus (+) icon will prompt you to view all the displays created by you or the dataset owner. To select/remove a specific display from this dataset, click the plus (+)/minus(-) button again. 

<img width="1374" height="648" alt="image" src="https://github.com/user-attachments/assets/e64e0bbb-118d-4d93-bacf-3008116b723e" />

---

**6.4 View and arrange the collection**

The layout of dataset collections (how they appear when viewed in the "Gene Expression Viewer") can be customized in the Dataset Explorer under "Collection arrangement view".

- In the "Dataset Explorer, select a collection: use the "Collection Management" dropdown to choose the dataset collection you want to edit. (Note: You can only modify collections you own. However, you may duplicate shared collections and then customize them as you wish). 
- Go to the "Collection arrangement view" by clicking the last icon among the view selection. If you do not have editing permissions for the collection, this button will not appear.
- Rearrange datasets by following the instructions. 
    - Drag and drop panels to change their order and position.
    - Resize panels by adjusting their width and height to make them larger, smaller, wider, or narrower.
    - When you are satisfied with the layout, click Save to apply your changes.

<img width="1374" height="648" alt="image" src="https://github.com/user-attachments/assets/8395bcc0-72e1-44b2-b113-e60471efbf5e" />

---

### 7. Further analysis using additional NeMO analytic features

**7.1 Gene List**

*7.1.1 Gene List Manager Overview*

The Gene List Manager in NeMO allows you to create, store, and organize user-defined sets of genes. The "Gene List Manager" follows a similar layout to the "Dataset Management", where you can search for other public or private gene lists, enabling easier access and comparison across multiple analyses. Instead of retyping gene names each time, you can save them into lists, quickly reuse them in different tools, and also share them with collaborators.

<img width="1374" height="648" alt="image" src="https://github.com/user-attachments/assets/f074ed15-aa64-4c4e-924e-474c84129777" />

1. Navigation & Access
- To begin creating a gene list, first open the gene list manager by clicking the "Gene List" tab in the navigation panel. 
- From here, you can access all your lists, lists shared with you, group-affiliated lists, or public lists.
- Lists can also be filtered by organism (e.g., human, mouse, zebrafish) and by date added.

2. Views
- At the top of the page, you can toggle between "List View", "Table View", and "Expanded View", similar to the "Dataset collection" a.
- These views control how gene lists are displayed, ranging from compact overviews to detailed, expanded records.
- Lists can be sorted by metadata such as date created. This is useful for finding your most recent lists quickly.

3. Core Functions for Each Gene List
Each entry in the "Gene list manager" represents a single gene list, showing its name, organism, owner, type, description, and added date. In the expanded gene list view entry, there are several action icons that allow you to manage the list:
- Projection Tool (eye icon) → Opens the list in NeMO’s projection tool for quick comparison across studies (see section 7.4 Projection Tool). 
- Download (down arrow icon) → Exports the list for offline use or backup.
- Delete (trash icon) → Permanently removes the list. 
- Shareable Link (link icon) → Generates a link to share the gene list with collaborators.
- Edit Link (pencil/edit icon) → Allows you to modify the sharing settings.
- Edit Metadata (page/edit icon) → Update metadata like the gene list’s description or category.

4. Clicking on the "+ Create new gene list" will lead you to the following page to set up your new gene list.  There are two types of gene lists NeMO users can create: unweighted and weighted. This tutorial 
There are two types of gene lists NeMO users can create: "unweighted gene list"(the most common type, just gene symbols), or "weighted gene list" (e.g., PCA loadings, differential gene statistics, etc.). We will cover how to create/upload them and use them, along with other NeMO analytic tools, in the next few sections.

<img width="1374" height="648" alt="image" src="https://github.com/user-attachments/assets/d00daad8-a2d0-4039-831d-d183055be73c" />

------

*7.1.2 Create unweighted gene list*

Following the instructions on the "uploading gene list, the users are allowed to create their gene list in two ways:

- Enter genes manually
<img width="1374" height="648" alt="image" src="https://github.com/user-attachments/assets/36c55d62-eb9c-4691-9ba4-b66839994c6f" />


- Upload a file following the format requirements
  - Select Upload file under Unweighted gene list.
  - Fill in the Name, Organism, and optional Description.
  - Click Choose a file… and upload a .txt, .csv, or .tsv file.
    - File must contain one gene symbol per line, or genes separated by commas, spaces, or tabs.
    - No headers or numeric weights should be included.

<img width="1374" height="612" alt="image" src="https://github.com/user-attachments/assets/ebb3956b-6e38-4597-8fad-b1fc009053bb" />

-----

*7.1.3 Create weighted gene list*

Similar to the unweighted gene list, users can also upload a weighted gene list in two different ways.
Weighted lists extend basic gene collections by associating each gene with one or more numerical weights — for example, PCA loadings, log₂ fold changes from differential expression, or gene enrichment scores.

This makes weighted gene lists particularly useful for:
- Visualizing ranked or gradient-based expression patterns
- Highlighting top contributing genes in dimensionality reduction (e.g., PCA, UMAP)
- Comparing differential analysis or enrichment results across datasets

Different from uploading the unweighted gene list file, the weighted gene list file has to have at least three columns
| Column                   | Description                                                                             |
| ------------------------ | --------------------------------------------------------------------------------------- |
| **1️ensembl_id**       | Unique gene identifier (required for cross-species mapping)                             |
| **2️gene_symbol**      | Human-readable gene name                                                                |
| **3️numeric columns** | One or more columns of numeric weights (e.g., PCA scores, DE log₂FC, enrichment values) |

<img width="1374" height="587" alt="image" src="https://github.com/user-attachments/assets/5eb86810-61a5-48f3-a86f-c30a1ac77e78" />


----
*7.1.4 Upload weighted gene list from Comparison tool (See 7.3.4)*

----

**7.2 Multi-gene viewer**

*7.2.1 Create multi-gene display*

Similar to the single-gene viewer, we can create different visualization panels for examining different groups of genes. The "multi-gene display Tool" allows users to visualize the expression of multiple genes simultaneously across defined biological or experimental conditions.
This function is ideal for comparing patterns of co-expression, examining pathway gene sets, or exploring marker profiles across disease, cell type, or genotype.

The multi-gene curator can be launched from the dataset management page under the analysis tool of each expanded dataset entry. Clicking on the "multi-gene curator" button will lead the user to the "multi-gene display page," and the selected dataset will be automatically populated in the "select a dataset" entry. Alternatively, the multi-gene curator can be launched using the navigation panel on the left. Using the drop-down and search button under the "select a dataset" entry, users can manually select their dataset of interest. 

<img width="1374" height="486" alt="image" src="https://github.com/user-attachments/assets/33fc760b-f4cb-423e-9503-577e6a758e47" />
<img width="1374" height="563" alt="image" src="https://github.com/user-attachments/assets/78713cb8-ce4c-4915-a903-f401da328cc6" />

Follow the prompts 2-5 in the " multi-gene display" page to customize the visualization for viewing multiple genes:
- Create or load an existing display
  - Click + Curate new display to start a new visualization, or select a saved display from your existing list
  - You can also Clone or Set as Default any previously created display for reuse
- Select plot type and analysis: these define how the gene expression data will be displayed and aggregated 
    - Add the genes you wish to visualize
    - Enter gene names manually into the search box, or use Quick search using Gene Lists to import a saved list.
    - As you type, the system will auto-suggest matching genes present in the dataset.
    - For multi-gene plots, at least two genes are required to create the plot template. Upon completion, the plot schematics can be applied to alternative genes or saved gene list when using the multi-gene display search tools. 
- Select plotting options: customize how your data will be grouped and displayed.
    - Primary grouping: choose the main variable for comparison (e.g., patho_diagnosis for AD vs. Control).
    - Secondary grouping (optional): add a secondary grouping factor if desired (e.g., cell type, genotype).
    - Add clusterbar groups: check boxes to add cluster bars (e.g., apoe_genotype, study, subset).
    - When all selections are complete, click "Plot" to generate your heatmap or matrix plot.

Clicking on the "plot" button will lead the user to the preview section, where plotting parameters can be further adjusted using the options listed on the left configuration panel. This page allows you to fine-tune the visualization and save it as a curated display. Once users are satisfied with the visualization, they can save the display and download it to their own device. 

<img width="1374" height="533" alt="image" src="https://github.com/user-attachments/assets/70dc2e1b-6992-45c4-8999-0ed0f20322d7" />

----

*7.2.2 Add the multi-gene display to a collection* 

As discussed in "6. Create dataset collections",  "single-gene display" can be added to a curated collection for cross-dataset comparison. Multi-gene displays can also be added and arranged in the existing collection or new collection in the same way. To do so, first navigate back to the "dataset management page". 

- Under Collection management, use the dropdown to:
    - Select an existing collection, or
    - Click + Create a new collection, then give it a name
    - Optionally, check the "Private collection" if you don’t want it shared publicly.

<img width="1374" height="711" alt="image" src="https://github.com/user-attachments/assets/e7277f4d-f2ee-4eec-b53f-b79463d02938" />

- Add the dataset containing your multi-gene display
    - Scroll through the dataset list (or use the keyword search/filter panel) to find the target dataset
    - Expand the dataset entry and use the "+/-" button to add or remove the display from the collection
<img width="1374" height="706" alt="image" src="https://github.com/user-attachments/assets/173b926f-989b-45d5-9e36-73630d030a3d" />


- Clicking on the "+/-" button will pop up a preview window to allow you to add or remove displays for this specific dataset
<img width="1374" height="326" alt="image" src="https://github.com/user-attachments/assets/6753c823-99bd-40a9-921c-3d5d0ce2cc9d" />

----

*7.2.3 View expression of genes saved in the weighted/unweighted genes*

Once displays are added to the collection, the multi-gene viewer function for this collection will be activated in the "Gene Expression" page. Go to the "Gene Expression" tab using the navigation panel on the left. 

- Select the collection that you added your multi-gene display to using the drop-down menu
- Type in the names of the genes or use the drop-down menu to select a saved gene list
  - Click to confirm the genes from the list
  - Confirm the query by clicking on proceed
- Run the query by clicking on the "magnifier" button at the top

<img width="1239" height="615" alt="image" src="https://github.com/user-attachments/assets/050e5897-d67c-4e59-afe5-fc7bb0359eb3" />

----

**7.3 Comparison tool**

The comparison tool is designed for comparing differences in gene expression between two groups of samples and testing for statistical differences in expression (Note: If the dataset does not have replicates, statistical tests will be unavailable.) It enables users to visualize differences in expression levels and fold-change, as well as to conduct statistical significance testing.

*7.3.1 Launch the comparison tool* 
There are two different ways to navigate to the comparison tool page:

1. To launch the comparison tool from the dataset collection page
   - Expand the dataset entry and click on the analysis tool. Use the drop-down menu to select the comparison tool. This will bring you to the comparison tool main page, and the selected dataset will be automatically populated on the comparison page.
   <img width="1374" height="565" alt="image" src="https://github.com/user-attachments/assets/00c236d5-8cf5-4fff-89ed-05f3545a0c73" />

2. Use the navigation panel and choose the dataset manually from the drop-down menu. Once selected, a green checkmark will confirm your dataset is loaded and ready for comparison.
   <img width="1374" height="533" alt="image" src="https://github.com/user-attachments/assets/fb6773a8-9dcb-4eb5-afab-e23935698da4" />

----

*7.3.2 Specify the  parameter*

- Condition parameter: 
    - Series to compare: choose the metadata field that represents your grouping variable (e.g., patho_diagnosis).
    - X-axis (query condition): select the condition to plot on the X-axis (e.g., AD).
    - Y-axis (reference condition): select the comparison or control group to plot on the Y-axis (e.g., Control).
    - You can include or exclude specific groups from the dropdown and toggle menus as needed. The number of observations (cells, nuclei, or samples) included in the comparison will appear below the selection fields.
<img width="1374" height="651" alt="image" src="https://github.com/user-attachments/assets/e6136dbf-7476-4123-b304-6b062c433050" />

 - Comparison parameter: to refine your analysis with statistical testing and filtering options. Click the Edit Parameters button at any time to revisit or adjust these settings.

    - Significance Test
      - Select test: choose a statistical test (e.g., t-test, Wilcoxon). If your dataset does not have replicates, this option will be greyed out (set to None).
      - P-value cutoff: specify the p-value threshold for significance (default: 0.05).
      - filter: decide whether to colorize significant points or filter out non-significant ones.
    
    - Data Filters
      - Report output as: Choose whether expression values should be reported as linear or log2.
      - Fold Change Cutoff (≥N): Set a threshold for the minimum fold change to display (default: 2.0). Genes below this cutoff are omitted from the plot for performance efficiency.
      - Deviation Filter: Optionally filter genes based on variability across samples (default: No filter).
<img width="1374" height="300" alt="image" src="https://github.com/user-attachments/assets/a61bf275-195a-488f-9151-5c924406b3e3" />

----

*7.3.3 Generate the Comparison Plot*

-  Once all the parameters are adjusted accordingly, select "Plot"
-  Each point represents a gene, plotted by its expression level in each condition:
    - X-axis: Expression in the query condition (e.g., AD)
    - Y-axis: Expression in the reference condition (e.g., Control)
- Interpretation:
    - Genes along the diagonal are similarly expressed in both conditions.
    - Genes deviating from the diagonal indicate differential expression.
    - Red points indicate genes passing the defined cutoffs (fold-change and/or significance), while gray points do not.
<img width="1374" height="572" alt="image" src="https://github.com/user-attachments/assets/1ebf6dcb-a717-4e40-ac29-6a06b2de0730" />

 ----
 
*7.3.4. Explore and Save Results* 

Once the plots are generated, users can interact with the plots for better visualization or save the results as a weighted gene list for projection or multi-gene display. 

- Use Plotly tools (upper right corner) to zoom, box-select, or lasso specific gene subsets.
    - Double-click to reset the selection.
    - Hover over points to view gene IDs and expression statistics.
<img width="1191" height="523" alt="image" src="https://github.com/user-attachments/assets/451b0f17-2c37-4290-bb8a-ac4bea11a32e" />

- Selected genes will appear in a table beneath the plot. You can save these into a gene list by entering a name and choosing:
    - Unweighted: saves only selected genes.
    - Weighted: saves all genes with their weights (log2 fold change values).
    - Click "Save" to store the list for later use in other NeMO analytical tools.
    - The saved gene list can be viewed from the "gene list manager" tab.
<img width="1135" height="667" alt="image" src="https://github.com/user-attachments/assets/a30f0a44-2529-48be-9f93-a507fbff8de4" />

- Highlight specific genes:
    - Search for specific genes using the search bar by typing the name of the gene. The selected genes will be highlighted on the plot or the table if there are selected clustered genes, using the "box select" or "lasso". 
    - Alternatively, use the Quick search using Gene Lists dropdown to highlight genes from existing gene lists.
<img width="1125" height="650" alt="image" src="https://github.com/user-attachments/assets/16b5a5a6-5754-488c-a9b2-31eb3197b891" />
<img width="1123" height="378" alt="image" src="https://github.com/user-attachments/assets/180812b8-4cf8-4da5-ad25-5ec35f68df43" />

----

**7.4 Projection tool**

*7.4.1 Overview of the projection tool*

High-dimensional genomics data often contain both biological signal and technical noise. Projection helps separate these components by transferring learned biological features from one dataset to another.
This allows researchers to validate biological patterns across independent datasets, identify shared or distinct expression programs across tissues, species, or disease states, and visualize how strongly each dataset expresses a given pattern or cell-state signature.
In short, projection provides a quantitative and visual way to interpret biological reproducibility and variation across studies. The Projection Tool in NeMO enables users to apply learned biological patterns from one dataset onto another, allowing comparison of gene expression dynamics across studies, cell types, or disease conditions. This approach builds on the principles of transfer learning introduced in the projectR framework (Stein-O’Brien et al., 2019), which uses dimension reduction outputs—such as PCA or NMF—to identify biologically meaningful axes that can be reused across related datasets.

*7.4.2 Projection in NeMO*

Within NeMO, projection is seamlessly integrated with your existing datasets and collections:

- Use the navigation panel on the left to launch the "Projection tool" 
- Search for a pattern
  - This pattern is derived from your "source" datasets.
  - This pattern could be saved unweighted gene list, a weighted gene list, or a weighted gene list saved from the comparison tool from your own collection of gene lists.
  - Alternatively, users can type the name of the gene list saved in the public domain of the NeMO gene list.
- Select a collection that contains visualization panels for your "target dataset" or where your pattern should be projected.

<img width="1374" height="668" alt="image" src="https://github.com/user-attachments/assets/8c4493b7-ada5-4011-939d-d41b4ff5ff27" />

- Select an algorithm — determine how the pattern is mathematically transferred
  - There are four different algorithms available on NeMO, which are shown in the table below.
  - In addition to the different algorithms, there are two options available for plotting the gene expression data.
    - Z-score normalize gene expression
      - When enabled, this option rescales each gene’s expression values across all samples so that the gene has a mean of 0 and standard deviation of 1. This normalization helps make genes more comparable by removing baseline expression level differences. However, it can also amplify noise
      - if the dataset contains many lowly expressed genes or dropout events, because small random variations get exaggerated after scaling.
    - Set negative coefficient weight outputs to zero
      - When enabled, any negative projection scores (e.g., negative pattern weights or component loadings) are set to zero during visualization.
      - This does not alter the saved projection data—it only affects what is shown in the plot.
      - This is useful for: making projected heatmaps or UMAPs easier to interpret when the focus is on positive pattern activation (e.g., enrichment of a biological program) and avoiding confusion from negative weights, which often represent anti-correlation rather than active expression.

| Algorithm                                   | Description                                                                                                          | Typical Use Case                                                        |
| ------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------- |
| **Principal Component Analysis (PCA)**      | Projects data onto principal axes, capturing the largest variance. Fast and generalizable.                            | Comparing overall gene expression trends between datasets.              |
| **Non-negative Matrix Factorization (NMF)** | Projects based on additive, parts-based gene programs. Useful for identifying cell-type or pathway-specific modules. | Cross-dataset validation of cell-type signatures or molecular programs. |
| **Fixed Gene Weights (from NMF re-run)**    | Applies existing NMF gene weights to a new dataset without recomputation.                                            | Efficient projection when using previously learned NMF signatures.      |
| **Binary Gene Count (planned)**             | Uses binary membership (gene present/absent) instead of continuous weights.                                          | Interpretable comparison of curated gene sets.                          |


- Use the "magnifier" button to confirm the projection parameter setting and start to run the projection
  - Each point in the resulting plot represents a cell or sample from the target dataset, colored by the strength of the projected signal.
  - If multiple patterns are present in the selected list, use different tabs to toggle between projections of different patterns.
  - Similar to the "single gene display" and the "multi-gene display" display. Users can enlarge the plots and switch between different displays from the same dataset using the "enlarge" button and the "for more" button.
<img width="1374" height="645" alt="image" src="https://github.com/user-attachments/assets/7baecdce-bbe7-4f5a-9654-ba1ac41051be" />
 
  - Clicking on "View all gene list" will open a separate window in which all the genes included in the pattern, with their corresponding loadings, are displayed
<img width="1374" height="645" alt="image" src="https://github.com/user-attachments/assets/3d579173-47e4-4637-8c1f-68703817c7c5" />
 

----

### 8. Reference
1. Shreyash Sonthalia, Ricky S. Adkins, Joshua Orvis, Guangyan Li, Xoel Mato Blanco, Alex Casella, Jinrui Liu, Genevieve Stein-O’Brien, Brian Caffo, Ronna Hertzano, Anup Mahurkar, Jesse Gillis, Jonathan Werner, Shaojie Ma, Nicola Micali, Nenad Sestan, Pasko Rakic, Gabriel Santpere, Seth A. Ament, Carlo Colantuoni.
A Curated Compendium of Transcriptomic Data for the Exploration of Neocortical Development.
bioRxiv (2024). doi:10.1101/2024.02.26.581612

2. Gaurav Sharma, Carlo Colantuoni, Loyal A. Goff, Elana J. Fertig, Genevieve Stein-O’Brien.
projectR: an R/Bioconductor package for transfer learning via PCA, NMF, correlation and clustering.
Bioinformatics (2020) 36(11): 3592–3593.

3. Gene Expression Analysis Resource (gEAR).
GitHub repository. University of Maryland School of Medicine.
https://github.com/IGS/gEAR
