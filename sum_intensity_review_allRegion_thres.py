import numpy as np
import tifffile
import os
import scipy
import dill
import pandas as pd

def load_object(filename):
    f = open(filename, "rb")
    obj = dill.load(f)
    f.close()
    return obj

## Directory of the folder containing activity images
file_dir="/home/sthyme/movedbutkeepinghere/ZZZSelectedImagespErk/autism_Oct2023/selectedmaintext_structure/"

## Load various masks
mask_noSC=load_object("/home/sthyme/movedbutkeepinghere/Zscripts_UAB_imageanalysis/whole_brain_mask_no_SpinalCord")
mask=load_object("/home/sthyme/movedbutkeepinghere/Zscripts_UAB_imageanalysis/whole_brain_mask")
all_masks_original=load_object("/home/sthyme/movedbutkeepinghere/Zscripts_UAB_imageanalysis/all_masks_sym_all_good")
gene_order=['arid1b','hectd4','kdm5b-ab','tlk2','phf2','nsd2-r378s','nsd2','deaf1-t238p','deaf1-c207y','deaf1','kiaa0232','cnot3-ab','tbl1xr1a','tcf7l2','phf12-ab','kmt2c-ab','gigyf1-ab','ash1l','kmt5b','hdlbp-ab']
selected_keys = ['Diencephalon___Dorsal_Thalamus','Diencephalon___Ventral_Thalamus','Diencephalon___Habenula','Diencephalon___Intermediate_Hypothalamus','Diencephalon___Rostral_Hypothalamus','Diencephalon___Caudal_Hypothalamus','Diencephalon___Posterior_Tuberculum','Diencephalon___Preoptic_Area','Diencephalon___Pretectum','Mesencephalon___Tectum_Neuropil','Mesencephalon___Tectum_Stratum_Periventriculare','Mesencephalon___Tegmentum','Mesencephalon___Torus_Longitudinalis','Mesencephalon___Torus_Semicircularis','Rhombencephalon___Cerebellum','Rhombencephalon','Telencephalon___Olfactory_Bulb','Telencephalon___Subpallium','Telencephalon___Pallium']


#dict_keys(['Diencephalon___Left_Habenula_Vglut2_Cluster', 'Mesencephalon___Torus_Longitudinalis', 'Telencephalon___Olfactory_Bulb', 'Rhombencephalon___Gad1b_Cluster_15', 'Rhombencephalon___Ventrolateral_population_of_serotonergic_neurons', 'Rhombencephalon___Medial_Vestibular_Nucleus', 'Rhombencephalon___Gad1b_Cluster_14', 'Rhombencephalon___Gad1b_Cluster_6', 'Mesencephalon___Otpb_Cluster', 'Rhombencephalon___Rhombomere_5', 'Rhombencephalon___Glyt2_Cluster_8', 'Telencephalon___Anterior_Commisure', 'Rhombencephalon___RoL2', 'Rhombencephalon___Isl1_Stripe_1', 'Diencephalon___Hypothalamus_Vglut2_Cluster_1', 'Rhombencephalon___Neuropil_Region_4', 'Rhombencephalon___Glyt2_Cluster_13', 'Rhombencephalon___RoL3', 'Rhombencephalon___67FDhcrtR_Gal4_Stripe_4', 'Diencephalon___Pineal_Vmat2_cluster', 'Mesencephalon___Vmat2_cluster2', 'Telencephalon___Isl1_cluster_1', 'Diencephalon___Retinal_Arborization_Field_6_AF6', 'Diencephalon___Preoptic_area_Vglut2_cluster', 'Rhombencephalon___Olig2_enriched_areas_in_cerebellum', 'Diencephalon___Dopaminergic_Cluster_3__hypothalamus', 'Rhombencephalon___Inferior_Olive', 'Rhombencephalon___Anterior_Cluster_of_nV_Trigeminal_Motorneurons', 'Rhombencephalon___Glyt2_Cluster_5', 'Rhombencephalon___Glyt2_Stripe_1', 'Telencephalon___Subpallial_Otpb_Cluster_2', 'Diencephalon___Caudal_Hypothalamus', 'Rhombencephalon___Qrfp_neuron_cluster_sparse', 'Rhombencephalon___Gad1b_Cluster_20', 'Ganglia___Lateral_Line_Neuromast_O1', 'Diencephalon___Oxtl_Cluster_4__sparse_in_hypothalamus', 'Diencephalon___Pretectal_dopaminergic_cluster', 'Rhombencephalon___Gad1b_Cluster_8', 'Mesencephalon___Tectum_Neuropil', 'Rhombencephalon___67FDhcrtR__Gal4_Stripe_1', 'Ganglia___Facial_glossopharyngeal_ganglion', 'Telencephalon___Vmat2_cluster', 'Rhombencephalon___Vglut2_Stripe_2', 'Diencephalon___Hypothalamus_Gad1b_Cluster_1', 'Rhombencephalon___Glyt2_Cluster_4', 'Mesencephalon___Isl1_cluster_of_the_mesencephalic_region', 'Rhombencephalon___Interpeduncular_Nucleus', 'Rhombencephalon___RoL__R1', 'Telencephalon___S1181t_Cluster', 'Diencephalon___Dopaminergic_Cluster_2__posterior_tuberculum', 'Diencephalon___Pretectum', 'Rhombencephalon___MiD2', 'Diencephalon___Hypothalamus_Olig2_cluster', 'Ganglia___Lateral_Line_Neuromast_D2', 'Diencephalon___Diffuse_Nucleus_of_the_Intermediate_Hypothalamus', 'Spinal_Cord___Isl1_stripe__motorneurons', 'Rhombencephalon___Glyt2_Cluster_10', 'Diencephalon___Anterior_group_of_the_posterior_tubercular_vmat2_neurons', 'Rhombencephalon', 'Mesencephalon___Vglut2_cluster_1', 'Diencephalon___Retinal_Arborization_Field_1_AF1', 'Rhombencephalon___Rhombomere_6', 'Rhombencephalon___X_Vagus_motorneuron_cluster', 'Rhombencephalon___Tangential_Vestibular_Nucleus', 'Rhombencephalon___Neuropil_Region_2', 'Mesencephalon___Retinal_Arborization_Field_7_AF7', 'Diencephalon___Dopaminergic_Cluster_4_5__posterior_tuberculum_and_hypothalamus', 'Rhombencephalon___Gad1b_Cluster_9', 'Spinal_Cord___Vmat2_Stripe_1', 'Rhombencephalon___Olig2_Stripe', 'Ganglia___Lateral_Line_Neuromast_OC1', 'Rhombencephalon___Gly2_Cluster_6', 'Diencephalon___Pretectal_Gad1b_Cluster', 'Diencephalon___Hypothalamus_Qrfp_neuron_cluster', 'Diencephalon___Preoptic_Area', 'Rhombencephalon___Mauthner', 'Rhombencephalon___MiV2', 'Diencephalon___Otpb_Cluster_2', 'Ganglia___Lateral_Line_Neuromast_D1', 'Diencephalon___Hypothalamus_Hcrt_Neurons', 'Rhombencephalon___Rhombomere_4', 'Diencephalon___Migrated_Area_of_the_Prectectum_M1', 'Rhombencephalon___Otpb_Cluster_3', 'Rhombencephalon___67FDhcrtR__Gal4_Cluster_4', 'Spinal_Cord___Vglut2_Stripe_3', 'Telencephalon___Pallium', 'Diencephalon___Anterior_preoptic_dopaminergic_cluster', 'Rhombencephalon___Gad1b_Cluster_4', 'Rhombencephalon___Glyt2_Stripe_2', 'Rhombencephalon___Vmat2_Cluster_4', 'Diencephalon___Preoptic_area_posterior_dopaminergic_cluster', 'Rhombencephalon___Glyt2_Cluster_9', 'Telencephalon___Postoptic_Commissure', 'Diencephalon___Optic_Chiasm', 'Diencephalon___Oxtl_Cluster_5', 'Rhombencephalon___S1181t_Cluster', 'Rhombencephalon___Eminentia_Granularis', 'Spinal_Cord___Vglut2_Stripe_1', 'Diencephalon___Right_Habenula_Vglu2_Cluster', 'Rhombencephalon___67FDhcrtR__Gal4_Cluster_2_Sparse', 'Rhombencephalon___Rhombomere_3', 'Rhombencephalon___RoM1', 'Rhombencephalon___Ptf1a_Cluster_1', 'Rhombencephalon___Otpb_Cluster_1', 'Rhombencephalon___Spiral_Fiber_Neuron_Anterior_cluster', 'Rhombencephalon___Cerebellum_Gad1b_Enriched_Areas', 'Diencephalon___Olig2_Band_2', 'Rhombencephalon___Isl1_Cluster_3', 'Diencephalon___Isl1_cluster_2', 'Rhombencephalon___Caudal_Ventral_Cluster_Labelled_by_Spinal_Backfills', 'Rhombencephalon___Mauthner_Cell_Axon_Cap', 'Ganglia___Lateral_Line_Neuromast_N', 'Ganglia___Trigeminal_Ganglion', 'Rhombencephalon___Vglut2_Stripe_3', 'Rhombencephalon___MiT', 'Diencephalon___Retinal_Arborization_Field_5_AF5', 'Diencephalon___Retinal_Arborization_Field_2_AF2', 'Telencephalon___Subpallium', 'Rhombencephalon___Rhombomere_7', 'Ganglia___Lateral_Line_Neuromast_SO3', 'Diencephalon___Otpb_Cluster_3', 'Diencephalon___Hypothalamic_VentroLateral_VMAT_cluster', 'Diencephalon___Oxtl_Cluster_1_in_Preoptic_Area', 'Rhombencephalon___Gad1b_Stripe_1', 'Rhombencephalon___Gad1b_Cluster_18', 'Rhombencephalon___Rhombomere_1', 'Diencephalon___Medial_vglut2_cluster', 'Telencephalon___Subpallial_Otpb_strip', 'Rhombencephalon___MiM1', 'Diencephalon___Posterior_Tuberculum', 'Rhombencephalon___Raphe__Superior', 'Mesencephalon___NucMLF_nucleus_of_the_medial_longitudinal_fascicle', 'Mesencephalon___Oxtl_Cluster_Sparse', 'Rhombencephalon___Glyt2_Stripe_3', 'Rhombencephalon___Gad1b_Cluster_3', 'Rhombencephalon___Raphe__Inferior', 'Rhombencephalon___Glyt2_Cluster_1', 'Rhombencephalon___Gad1b_Cluster_19', 'Diencephalon___Hypothalamus_6.7FRhcrtRGal4_cluster_1', 'Mesencephalon___Torus_Semicircularis', 'Rhombencephalon___Lobus_caudalis_cerebelli', 'Rhombencephalon___MiR2', 'Diencephalon___Hypothalamus_Vglut2_Cluster_5', 'Rhombencephalon___Vglut2_cluster_2', 'Rhombencephalon___RoM3', 'Rhombencephalon___MiD3', 'Rhombencephalon___Oculomotor_Nucleus_nIV', 'Rhombencephalon___Valvula_Cerebelli', 'Spinal_Cord___67FDhcrtR__Gal4_Stripe', 'Rhombencephalon___Cerebellum', 'Mesencephalon___Vmat2_cluster_of_paraventricular_organ', 'Rhombencephalon___67FDhcrtR__Gal4_Cluster_5', 'Diencephalon___Oxtl_Cluster_2', 'Spinal_Cord___Neuropil_Region', 'Diencephalon___Hypothalamus_Vglut2_Cluster_2', 'Diencephalon___Preoptic_Otpb_Cluster', 'Rhombencephalon___Gad1b_Cluster_12', 'Rhombencephalon___Vmat2_Cluster_3', 'Rhombencephalon___Small_cluster_of_TH_stained_neurons', 'Rhombencephalon___Vmat2_Cluster_1', 'Diencephalon___Retinal_Arborization_Field_3_AF3', 'Rhombencephalon___Rhombomere_2', 'Telencephalon___Telencephalic_Migrated_Area_4_M4', 'Telencephalon___Vglut2_rind', 'Rhombencephalon___Gad1b_Stripe_3', 'Rhombencephalon___Glyt2_Cluster_11', 'Rhombencephalon___Spiral_Fiber_Neuron_Posterior_cluster', 'Rhombencephalon___Locus_Coreuleus', 'Rhombencephalon___Gad1b_Cluster_10', 'Rhombencephalon___Oxtl_Cluster_2_near_MC_axon_cap', 'Mesencephalon___Oculomotor_Nucleus_nIII', 'Ganglia___Lateral_Line_Neuromast_SO1', 'Diencephalon', 'Diencephalon___Hypothalamus_Olig2_cluster_2', 'Mesencephalon___Ptf1a_Cluster', 'Rhombencephalon___Neuropil_Region_5', 'Ganglia___Eyes', 'Diencephalon___Torus_Lateralis', 'Diencephalon___Retinal_Arborization_Field_4_AF4', 'Rhombencephalon___Glyt2_Cluster_7', 'Diencephalon___Habenula', 'Spinal_Cord___Gad1b_Stripe_1', 'Rhombencephalon___Glyt2_Cluster_12', 'Rhombencephalon___67FDhcrtR_Gal4_Stripe_3', 'Rhombencephalon___VII_Facial_Motor_and_octavolateralis_efferent_neurons', 'Rhombencephalon___Gad1b_Cluster_7', 'Mesencephalon___Sparse_67FRhcrtR_cluster', 'Rhombencephalon___Posterior_Cluster_of_nV_Trigeminal_Motorneurons', 'Diencephalon___Oxtl_Cluster_3', 'Diencephalon___Hypothalamus__Caudal_Hypothalamus_Neural_Cluster', 'Rhombencephalon___Vmat2_Stripe_1', 'Diencephalon___Pineal', 'Telencephalon___Olfactory_bulb_dopaminergic_neuron_areas', 'Diencephalon___Hypothalamus_6.7FRhcrtRGal4_cluster_2', 'Rhombencephalon___Vmat2_Cluster_5', 'Diencephalon___Hypothalamus_s1181t_Cluster', 'Rhombencephalon___Corpus_Cerebelli', 'Mesencephalon___Retinal_Arborization_Field_8_AF8', 'Rhombencephalon___Gad1b_Cluster_11', 'Diencephalon___Otpb_Cluster_1', 'Rhombencephalon___Neuropil_Region_3', 'Diencephalon___Isl1_cluster_3', 'Rhombencephalon___Gad1b_Cluster_17', 'Rhombencephalon___RoV3', 'Ganglia___Lateral_Line_Neuromast_SO2', 'Rhombencephalon___67FDhcrtR__Gal4_Cluster_3', 'Rhombencephalon___Gad1b_Cluster_16', 'Mesencephalon', 'Rhombencephalon___Otpb_Cluster_4', 'Rhombencephalon___Lateral_Reticular_Nucleus', 'Ganglia___Olfactory_Epithelium', 'Diencephalon___Olig2_Band', 'Diencephalon___Otpb_Cluster_4', 'Telencephalon___Subpallial_dopaminergic_cluster', 'Spinal_Cord___Gad1b_Stripe_2', 'Diencephalon___Ventral_Thalamus', 'Spinal_Cord', 'Rhombencephalon___Vmat2_Cluster_2', 'Telencephalon___Subpallial_Gad1b_Cluster', 'Rhombencephalon___Oxtl_Cluster_1_Sparse', 'Rhombencephalon___CaV', 'Rhombencephalon___Area_Postrema', 'Spinal_Cord___Vglut2_Stripe_2', 'Rhombencephalon___Cerebellar__Vglut2_enriched_areas', 'Rhombencephalon___Spinal_Backfill_Vestibular_Population', 'Spinal_Cord___Olig2_Stripe', 'Diencephalon___Dopaminergic_Cluster_7__Caudal_Hypothalamus', 'Telencephalon___Subpallila_Vglut2_Cluster', 'Rhombencephalon___Vmat2_Stripe_3', 'Rhombencephalon___67FDhcrtR__Gal4_Cluster_1', 'Diencephalon___Dorsal_Thalamus', 'Rhombencephalon___Cerebellar_Neuropil_1', 'Spinal_Cord___Glyt2_Stripe', 'Rhombencephalon___67FDhcrtR__Gal4_Stripe_2', 'Rhombencephalon___Gad1b_Stripe_2', 'Ganglia___Posterior_Lateral_Line_Ganglia', 'Rhombencephalon___Gad1b_Cluster_2', 'Spinal_Cord___Vmat2_Stripe_2', 'Diencephalon___Dopaminergic_Cluster_6__hypothalamus', 'Rhombencephalon___Glyt2_Cluster_2', 'Telencephalon___Olig2_Cluster', 'Rhombencephalon___Vglut2_cluster_1', 'Diencephalon___Pituitary', 'Spinal_Cord___Dorsal_Sparse_Isl1_cluster', 'Diencephalon___Isl1_cluster_1', 'Diencephalon___Intermediate_Hypothalamus', 'Diencephalon___Hypothalamus_Vglut2_Cluster_6', 'Rhombencephalon___Gad1b_Cluster_13', 'Ganglia___Vagal_Ganglia', 'Ganglia___Facial_Sensory_Ganglion', 'Mesencephalon___Tectum_Stratum_Periventriculare', 'Rhombencephalon___VII_prime_Facial_Motor_and_octavolateralis_efferent_neurons', 'Rhombencephalon___Neuropil_Region_6', 'Mesencephalon___Retinal_Arborization_Field_9_AF9', 'Rhombencephalon___RoM2', 'Rhombencephalon___Isl1_Cluster_1', 'Diencephalon___Migrated_Posterior_Tubercular_Area_M2', 'Rhombencephalon___CaD', 'Telencephalon___Isl1_cluster_2', 'Rhombencephalon___MiR1', 'Ganglia___Anterior_Lateral_Line_Ganglion', 'Diencephalon___Hypothalamus_Gad1b_Cluster_3_Sparse', 'Rhombencephalon___Vmat2_Stripe_2', 'Telencephalon', 'Diencephalon___Hypothalamus_Vglut2_Cluster_3', 'Rhombencephalon___Glyt2_Cluster_3', 'Rhombencephalon___Olig2_Cluster', 'Telencephalon___Optic_Commissure', 'Rhombencephalon___Glyt2_Cluster_14', 'Diencephalon___Dopaminergic_Cluster_1__ventral_thalamic_and_periventricular_posterior_tubercular_DA_neurons', 'Rhombencephalon___MiV1', 'Rhombencephalon___Vglut2_Stripe_1', 'Rhombencephalon___Vglut2_Stripe_4', 'Ganglia___Statoacoustic_Ganglion', 'Rhombencephalon___Otpb_Cluster_5', 'Rhombencephalon___Gad1b_Cluster_5', 'Diencephalon___Hypothalamus_Gad1b_Cluster_2', 'Mesencephalon___Tegmentum', 'Rhombencephalon___Otpb_Cluster_6', 'Diencephalon___Anterior_pretectum_cluster_of_vmat2_neurons', 'Rhombencephalon___Gad1b_Cluster_1', 'Diencephalon___Rostral_Hypothalamus', 'Diencephalon___Hypothalamus_Vglut2_Cluster_4', 'Diencephalon___Eminentia_Thalami', 'Rhombencephalon___Ptf1a_Stripe', 'Spinal_Cord___Neurons_with_descending_projections_labelled_by_spinal_backfills', 'Rhombencephalon___Isl1_Cluster_2', 'Rhombencephalon___Vglut2_cluster_3', 'Diencephalon___Hypothalamus__Interediate_Hypothalamus_Neural_Cluster', 'Rhombencephalon___Noradrendergic_neurons_of_the_interfascicular_and_Vagal_areas', 'Rhombencephalon___Vglut2_cluster_4', 'Mesencephalon___Medial_Tectal_Band', 'Rhombencephalon___Otpb_Cluster_2__locus_coeruleus', 'Diencephalon___Postoptic_Commissure'])


all_masks = {k: v for k, v in all_masks_original.items() if k in selected_keys}

#print(all_masks_original.keys())
#exit()
#mask_dien=all_masks['Diencephalon']
#mask_dthal=all_masks['Diencephalon__Dorsal_Thalamus']
#mask_vthal=all_masks['Diencephalon__Ventral_Thalamus']
#mask_hab=all_masks['Diencephalon__Habenula']
#mask_ihyp=all_masks['Diencephalon___Intermediate_Hypothalamus']
#mask_pt=all_masks['Diencephalon___Posterior_Tuberculum']
#mask_mesen=all_masks['Mesencephalon']
#mask_rhom=all_masks['Rhombencephalon']
#mask_telen=all_masks['Telencephalon']
#mask_pal=all_masks['Telencephalon__pallium']
#mask_spal=all_masks['Telencephalon__subpallium']
#print(all_masks.keys())
#exit()

## Parse image names to extract genotype information: gene name, hom/het
files=[name for name in os.listdir(file_dir)]
genes=np.array([name.split("_")[0] for name in files])
#genos=np.array([name.split("_")[2] for name in files])
#genos=np.array([''.join(i for i in j if not i.isdigit()) for j in genos])
labels=[genes[i] for i in np.arange(len(genes))]
#labels=[genes[i]+"_"+genos[i] for i in np.arange(len(genes))]

## Preallocate dataframes (1 for each channel) to store sum of intensities in each image (rows) for each region (columns)
all_regions=list(all_masks.keys())
sum_G=pd.DataFrame(np.zeros((len(files)+1,len(all_regions))),index=["size"]+labels,columns=all_regions)
sum_R=pd.DataFrame(np.zeros((len(files)+1,len(all_regions))),index=["size"]+labels,columns=all_regions)
#sum_both=pd.DataFrame(np.zeros((len(files)+1,len(all_regions))),index=["size"]+labels,columns=all_regions)

## Calculate the size of each brain region by summing up each region mask. Write the sums in the dataframes as a row
sum_G.loc["size",:]=[np.sum(all_masks[region_mask]) for region_mask in all_regions]
sum_R.loc["size",:]=[np.sum(all_masks[region_mask]) for region_mask in all_regions]
#sum_both.loc["size",:]=[np.sum(all_masks[region_mask]) for region_mask in all_regions]

## Preallocate dataframes to store number of active pixels in each image for each region
sum_numG=pd.DataFrame(np.zeros((len(files)+1,len(all_regions))),index=["size"]+labels,columns=all_regions)
sum_numR=pd.DataFrame(np.zeros((len(files)+1,len(all_regions))),index=["size"]+labels,columns=all_regions)
#sum_numboth=pd.DataFrame(np.zeros((len(files)+1,len(all_regions))),index=["size"]+labels,columns=all_regions)

sum_numG.loc["size",:]=[np.sum(all_masks[region_mask]) for region_mask in all_regions]
sum_numR.loc["size",:]=[np.sum(all_masks[region_mask]) for region_mask in all_regions]
#sum_numboth.loc["size",:]=[np.sum(all_masks[region_mask]) for region_mask in all_regions]

#labels=[filename.split('_')[0] for filename in files]

#sum_G=pd.DataFrame(np.zeros((len(files),6)),index=labels,columns=['Brain','NoSpinalCord','Diencephalon','Mesencephalon','Rhombencephalon','Telencephalon'])
#sum_R=pd.DataFrame(np.zeros((len(files),6)),index=labels,columns=['Brain','NoSpinalCord','Diencephalon','Mesencephalon','Rhombencephalon','Telencephalon'])
#sum_both=pd.DataFrame(np.zeros((len(files),6)),index=labels,columns=['Brain','NoSpinalCord','Diencephalon','Mesencephalon','Rhombencephalon','Telencephalon'])

## set intensity threshold for calling active pixels. 
thres=10
## Calculate region-wise sum of intensities and number of active pixels for each image
for i in np.arange(len(files)):
    file_name=files[i]
    label=labels[i]
    print("summing intensities for "+label+"...")
    brain_im=np.array(tifffile.imread(file_dir+file_name))
    brain_R=brain_im[:,:,:,0]
    brain_R=brain_R*(brain_R>=thres)
    brain_G=brain_im[:,:,:,1]
    brain_G=brain_G*(brain_G>=thres)
    #brain_both=np.max(np.array([brain_im[:,:,:,0],brain_im[:,:,:,1]]),axis=0)
    #sum_G.loc[label,:]=[np.sum(brain_G*mask),np.sum(brain_G*mask_noSC),np.sum(brain_G*mask_dien),np.sum(brain_G*mask_mesen),np.sum(brain_G*mask_rhom),np.sum(brain_G*mask_telen)]
    #sum_R.loc[label,:]=[np.sum(brain_R*mask),np.sum(brain_R*mask_noSC),np.sum(brain_R*mask_dien),np.sum(brain_R*mask_mesen),np.sum(brain_R*mask_rhom),np.sum(brain_R*mask_telen)]
    #sum_both.loc[label,:]=[np.sum(brain_both*mask),np.sum(brain_both*mask_noSC),np.sum(brain_both*mask_dien),np.sum(brain_both*mask_mesen),np.sum(brain_both*mask_rhom),np.sum(brain_both*mask_telen)]
    sum_G.loc[label,:]=[np.sum(brain_G*all_masks[region_mask]) for region_mask in all_regions]
    sum_R.loc[label,:]=[np.sum(brain_R*all_masks[region_mask]) for region_mask in all_regions]
    #sum_both.loc[label,:]=[np.sum(brain_both*all_masks[region_mask]) for region_mask in all_regions]
    
    sum_numG.loc[label,:]=[np.sum((brain_G>0)*all_masks[region_mask]) for region_mask in all_regions]
    #sum_numG = sum_numG.add_suffix('_aaaGreen')
    sum_numR.loc[label,:]=[(np.sum((brain_R>0)*all_masks[region_mask])) for region_mask in all_regions]
    #sum_numR.loc[label,:]=[-1*(np.sum((brain_R>0)*all_masks[region_mask])) for region_mask in all_regions]
    #sum_numR = sum_numR.add_suffix('_aaaRed')
    #sum_numboth.loc[label,:]=[np.sum((brain_both>0)*all_masks[region_mask]) for region_mask in all_regions]


## save the dataframes
#sum_G.to_csv('./all_regions_sum_perk_green_channel_PaperData_thres50.csv')
#sum_R.to_csv('./all_regions_sum_perk_red_channel_PaperData_thres50.csv')
#sum_both.to_csv('/n/schier_lab2/users/yiqunwang/Summer Data/ReviewAnalysis/intensity_sum/all_regions_sum_perk_both_channels_PaperData_thres50.csv')

size_row = sum_numG.loc["size"]
sum_numG = sum_numG.div(size_row, axis=1)
sum_numR = sum_numR.div(size_row, axis=1)
#print(sum_numG)
sum_numG = sum_numG.reindex(columns=selected_keys)
sum_numR = sum_numR.reindex(columns=selected_keys)
sum_numG = sum_numG.reindex(index=gene_order)
sum_numR = sum_numR.reindex(index=gene_order)

#sum_numG = sum_numG[selected_keys].sort_index(axis=1)
#print(sum_numG)
#sum_numR = sum_numR[selected_keys].sort_index(axis=1)
sum_numG.to_csv('./structurefixjan2024_all_regions_sum_nPix_perk_green_channel_PaperData_thres50.csv')
sum_numR.to_csv('./structurefixjan2024_all_regions_sum_nPix_perk_red_channel_PaperData_thres50.csv')
#sum_numboth.to_csv('/n/schier_lab2/users/yiqunwang/Summer Data/ReviewAnalysis/intensity_sum/all_regions_sum_nPix_perk_both_channels_PaperData_thres50.csv')

#comboDF = pd.concat([sum_numG, sum_numR], axis=1)
#comboDF = comboDF.sort_index(axis=1)
#comboDF.to_csv('./novredo_all_regions_sum_nPix_perk_two_channel_PaperData_thres50.csv')
