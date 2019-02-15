# -*- coding: UTF-8 -*-
# create_paper_figures.py

import numpy as np
import pandas as pd
import statsmodels.formula.api as sm
import os

import cantera_tools as ctt
import analysis_methods as am

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib.lines import Line2D

image_path = '../results'
if not os.path.exists(image_path):
    os.makedirs(image_path)

# set plot style
sns.set_palette('colorblind',n_colors=4)
sns.set_style('white')
sns.set_context('paper',font_scale=1.5)
sns.set_style('ticks',{'ytick.direction': 'in','xtick.direction': 'in'})
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

# get data for models
# cluster number
model_molecule_to_cluster_number = {'full': {'propane':26, 'ethyl':25,'methyl':24,
                                             'ethene':22,'H-atom':23,'n-propyl':20,
                                             'methane':19,'ethane':18,'ethenyl':17,
                                             'ethyne':16,'propene':14,},
                                    'drg' : {'propane':8, 'ethyl':6,'methyl':5,
                                             'ethene':2,'H-atom':4,'n-propyl':3,
                                             'methane':1,'ethane':7,},
                                    '3rxn' : {'propane':5, 'ethyl':3,'methyl':2,
                                             'ethene':1,
                                             'methane':0,'ethane':4,},
                                    '6rxn' : {'propane':7, 'ethyl':5,'methyl':4,
                                             'ethene':2,'H-atom':3,
                                             'methane':1,'ethane':6,},
                                    }
# initial propane isotopologue concentrations
delta_total= -28
psia = 5.4
edge_labeled_delta = delta_total + psia / 3.
center_labeled_delta = delta_total - 2. * psia / 3.
edge_labeled_fraction = am.getEnrichementFractionFromDelta(edge_labeled_delta)
center_labeled_fraction = am.getEnrichementFractionFromDelta(center_labeled_delta)

fraction_propane = 0.0049 # see supplemental

initialMoleFractions={
"CCC": fraction_propane * (1-center_labeled_fraction) * (1-edge_labeled_fraction)**2,
"CCC-2": fraction_propane * center_labeled_fraction * edge_labeled_fraction**2,
"CCC-3": fraction_propane * edge_labeled_fraction**2*(1-center_labeled_fraction),
"CCC-4": fraction_propane * 2*edge_labeled_fraction * (1-edge_labeled_fraction) * center_labeled_fraction,
"CCC-5": fraction_propane * 2*edge_labeled_fraction *(1-center_labeled_fraction) * (1-edge_labeled_fraction),
"CCC-6": fraction_propane * center_labeled_fraction*(1-edge_labeled_fraction)**2,
"[He]": 1-fraction_propane,
}
main_paths = [('full', '../mechanisms/full_model'),
              ('drg','../mechanisms/drg_model'),
              ('3rxn','../mechanisms/three_reaction_model'),
              ('6rxn','../mechanisms/six_reaction_model')]

################################
# Figures 2 and 3
print('creating figures 2 and 3')

enrichment_results = []
concentrations = {}
ethyl_psie_all = pd.DataFrame()

# run all four simulations
for name, mainPath in main_paths:
    cluster_info = pd.read_csv(os.path.join(mainPath, 'isotopomer_cluster_info.csv'),index_col='name')
    
    molecule_to_cluster_number = model_molecule_to_cluster_number[name]

    temp = 850+273
    times = np.linspace(1e-4,95. / temp,100)

    solution = ctt.create_mechanism(os.path.join(mainPath,'chem.cti'))

    conditions = temp, 2e5, initialMoleFractions

    output = ctt.run_simulation(solution, times, conditions=conditions, 
                          condition_type = 'constant-temperature-and-pressure',
                          output_species = True,
                          output_reactions = False,
                          output_directional_reactions = False,
                          output_rop_roc = False)
    species = output['species']

    # find enrichments and total concentration
    delta_enrichments = pd.DataFrame(columns=list(molecule_to_cluster_number.keys()), index = times)
    concentration_data = pd.DataFrame(columns=list(molecule_to_cluster_number.keys()), index = times)
    
    ethyl_psie = pd.Series()
    for time in times:
        enrichments_temp = {}
        for molecule, cluster_num in molecule_to_cluster_number.items():
            labels = cluster_info[cluster_info.cluster_number == cluster_num].index
            if not np.isclose(species.loc[time,labels].sum(),0,atol=1e-40) and molecule != 'H-atom':
                enrichments_temp[molecule] = am.getDelta(species.loc[time,labels],
                                           cluster_info.loc[labels,'enriched_atoms'],
                                           cluster_info.loc[labels,'unenriched_atoms'],)
            concentration_data.loc[time,molecule] = species.loc[time,labels].sum()
            # get psie information
            if molecule == 'ethyl':
                ethyl_psie[time] = am.getPSIE(species.loc[time,labels],
                                               cluster_info.loc[labels,:],
                                               type_1 = 'r',
                                               type_2 = 'not_r')
        delta_enrichments.loc[time,:] = enrichments_temp
    concentrations[name] = concentration_data
    ethyl_psie_all[name] = ethyl_psie
    enrichment_results.append((name,delta_enrichments))

# output figure 2

molecule = 'ethene'

ax = plt.subplot()
for identifier, enrichments in enrichment_results:
    enrichments.plot(y=molecule,ax=ax)
    
ax.set_ylabel(u'{}\n$^{{13}}\delta$C\n(‰)'.format(molecule),rotation='horizontal',labelpad=30)
ax.set_xlabel('time (s)')

ax.annotate('6 reaction model',(.22,.4),xycoords='axes fraction', rotation = 0)
ax.annotate('full model',(.7,.7),xycoords='axes fraction', rotation = 3.5)
ax.annotate('3 reaction model',(.2,.88),xycoords='axes fraction', rotation = 3.5)
ax.annotate('DRG model',(0.7,.82),xycoords='axes fraction', rotation = 3)

ax.legend([])
plt.savefig(os.path.join(image_path,'{}_enrich.pdf'.format(molecule)),bbox_inches='tight')

# output figure 3

molecule = 'ethyl'

ax = plt.subplot()
for column_name in ethyl_psie_all.columns:
    psie_data = ethyl_psie_all[column_name]
    psie_data.plot(ax=ax)
    
ax.set_ylabel(u'{}\nPSIA\n(‰)'.format(molecule),rotation='horizontal',labelpad=20)
ax.set_xlabel('time (s)')

ax.annotate('6 reaction model',(.15,.95), (.2,.7),xycoords='axes fraction',textcoords='axes fraction', rotation = 0, arrowprops={'arrowstyle':'-'}) #kie
ax.annotate('full model',(.4,.91), (.4,.8),xycoords='axes fraction',textcoords='axes fraction', rotation = 0, arrowprops={'arrowstyle':'-'})
ax.annotate('DRG model',(.65,.89), (.6,.7),xycoords='axes fraction',textcoords='axes fraction', rotation = 0, arrowprops={'arrowstyle':'-'})
ax.annotate('3 reaction model',(.2,.5),xycoords='axes fraction', rotation = 0)

ax.set_ylim(-21,1)
ax.legend([])
plt.savefig(os.path.join(image_path,'{}_psie.pdf'.format(molecule)),bbox_inches='tight')

##############################################
# figure 5 - mole fractions
print('creating figure 4')
# get data
model_index = 0
model_name = main_paths[model_index][0]
model_path = main_paths[model_index][1]
molecule_to_cluster_number = model_molecule_to_cluster_number[model_name]
model_file = os.path.join(model_path, 'chem.cti')
isotopomer_info = pd.read_csv(os.path.join(model_path,'isotopomer_cluster_info.csv'),index_col = 'name')

# simulate mechanism
mole_fractions = pd.DataFrame(index=molecule_to_cluster_number.keys())
for temp in (750,800,850,900,950):

    solution = ctt.create_mechanism(model_file)

 
    # set initial conditions of solution in kelvin pascals and mole fractions
    conditions = temp+273, 2e5, initialMoleFractions
    t_final = 95. / temp # see supplemental info

    output = ctt.run_simulation(solution, [0,t_final], conditions=conditions, 
                          condition_type = 'constant-temperature-and-pressure',
                          output_species = True,
                          output_reactions = True,
                          output_directional_reactions = False,
                          output_rop_roc = False)
    species = output['species'].loc[t_final,:]

    isotopomer_info['conc'] = species

    # Gilbert et al 2016 weighted the results by the area of detector peak, 
    # so the total concentration should be weighted
    # by the number of carbon atoms.
    species_concentrations = {}
    for molecule, cluster_num in molecule_to_cluster_number.items():
        labels = isotopomer_info[isotopomer_info.cluster_number == cluster_num].index
        species_concentrations[molecule] = isotopomer_info.loc[labels,'conc'].sum() * \
                                                (isotopomer_info.loc[labels[0],'enriched_atoms'] + isotopomer_info.loc[labels[0],'unenriched_atoms'])

    # find mole fractions (weighted by number of carbons)
    total_concentration = 0 
    for index in isotopomer_info.index:
        moles_element = isotopomer_info.loc[index,'conc'] * (isotopomer_info.loc[index,'enriched_atoms'] + isotopomer_info.loc[index,'unenriched_atoms'])
        total_concentration += moles_element

    mole_fractions_temp = {}  
    for molecule in species_concentrations.keys():
        mole_fractions_temp[molecule] = species_concentrations[molecule] / total_concentration
    mole_fractions[temp] = pd.Series(mole_fractions_temp)

# get experimental data
# taken three times to improve accuracy
experimental_data_folder = '../exp_data/'
fig1A_data_a = pd.read_csv(os.path.join(experimental_data_folder,'Gilbert_fig1A_from_engauge_try_2a.csv'), 
                           index_col='Temperature (C)')
fig1A_data_b = pd.read_csv(os.path.join(experimental_data_folder,'Gilbert_fig1A_from_engauge_try_2b.csv'), 
                           index_col='Temperature (C)')
fig1A_data_c = pd.read_csv(os.path.join(experimental_data_folder,'Gilbert_fig1A_from_engauge_try_2c.csv'), 
                           index_col='Temperature (C)')
fig1A_data_original = (fig1A_data_a + fig1A_data_b + fig1A_data_c)/3

# process  data
mole_fractions_1a = mole_fractions.T * 100
# normalize using propane conversion for model. use 100% for experiments
mole_fractions_1a = mole_fractions_1a.divide(100-mole_fractions_1a['propane'],'index') * 100
column_order = ['methane', 'ethane','ethene','propene']
fig1A_data = fig1A_data_original[column_order]
column_order = ['methane', 'ethane','ethene','propene','ethyne']
mole_fractions_1a = mole_fractions_1a[column_order]
fig1A_data = fig1A_data.divide(fig1A_data.sum('columns'),'index')*100

# create figure 4

f,ax = plt.subplots()

fig1A_data.plot.area(ax=ax, linewidth=0, 
                stacked=True, alpha= 0.3)
mole_fractions_1a.plot(ax=ax, linestyle = '-', 
                       linewidth= 2, markersize =0,stacked=True)
ax.set_ylabel("fraction carbon from propane (%)")
ax.set_xlabel(u'T (°C)')
ax.set_xticks(ticks=[750,800,850,900,950])
ax.set_xlim(750,950)
ax.set_ylim(1,110)

ax.legend([])
ax.annotate('methane',(800,18),(760,5),arrowprops={'arrowstyle':'-'})
ax.annotate('ethane',(800,31),(760,22),arrowprops={'arrowstyle':'-'})
ax.annotate('ethene',(800,91),(760,55),arrowprops={'arrowstyle':'-'})
ax.annotate('propene',(820,97),(760,92),arrowprops={'arrowstyle':'-'})
ax.annotate('ethyne',(900,100),(860,103),arrowprops={'arrowstyle':'-'})

plt.savefig(os.path.join(image_path, '1a_area_normalized_using_experimental_conversion.pdf'), bbox_inches = 'tight')

################################
# figure 5 - enrichments
print('creating figure 5')
# plot experimental data
# simulate and plot model data
styles = ['o',
         '^',
         's',
         'v',
         'D']
line_styles= [(0, (1, 2)),
             (0, (5, 10)),
             (0, (1, 5)),
             (0, (3, 5, 1, 5)),
             ]
f,ax = plt.subplots()
fig1B_data = pd.read_csv('../exp_data/Gilbert_fig1B_engauge.csv', index_col='Temperature (C)')
del fig1B_data['propene']
fig1B_data.plot(ax=ax, linestyle = '',
                linewidth= .5, style = styles, markersize = 6, markeredgewidth = 1)



# use same enrichments as used by gilbert
conversions_gilbert = 1- fig1A_data_original.propane / 100

for model_index in range(2):
    # get model data
    model_name = main_paths[model_index][0]
    model_path = main_paths[model_index][1]
    molecule_to_cluster_number = model_molecule_to_cluster_number[model_name]
    model_file = os.path.join(model_path, 'chem.cti')
    isotopomer_info = pd.read_csv(os.path.join(model_path,'isotopomer_cluster_info.csv'),index_col = 'name')

    # simulate model
    enrichments_by_conversion = pd.DataFrame(dtype=float)
    for temp in [750,800,850,900,950]:
        conversion = conversions_gilbert[temp]

        # creates the cantera Solution object
        solution = ctt.create_mechanism(model_file)

        conditions = temp+273, 2e5, initialMoleFractions

        output = ctt.run_simulation_till_conversion(solution, species='CCC', conversion=conversion,conditions = conditions,
                              condition_type = 'constant-temperature-and-pressure',
                              output_species = True,
                              output_reactions = False,
                              output_directional_reactions = False,
                              output_rop_roc = False)
        species = output['species'].iloc[-1,:]

        isotopomer_info['conc'] = species

        for molecule, cluster_num in molecule_to_cluster_number.items():
            labels = isotopomer_info[isotopomer_info.cluster_number == cluster_num].index
            if molecule != 'H-atom':
                enrichments_by_conversion.loc[temp, molecule] = am.getDelta(species[labels],
                                                   isotopomer_info.loc[labels,'enriched_atoms'],
                                                   isotopomer_info.loc[labels,'unenriched_atoms'],
                                                   )
    # plot this data set
    enrichments_by_conversion.plot(y=fig1B_data.columns,
                  ax=ax, linestyle = line_styles[model_index], linewidth= 1, 
                  style =styles,markersize = 2, markeredgewidth = 1, 
                  markerfacecolor = "None")
# plot details
ax.set_ylabel("$\delta^{13}C$\n$(\perthousand)$",rotation='horizontal',va='center',ha='right')
ax.set_xlabel(u'T (°C)')
# move legend outside of plot
ax.legend(list(fig1B_data.columns), **{'bbox_to_anchor':(1.05, 1), 'loc':2, 'borderaxespad':0.})

items, entries = ax.get_legend_handles_labels()
items = items[:len(fig1B_data.columns)]
legend_items = []
for item in items:
    # make a new copy so the graph isn't affected without the lines
    legend_item = Line2D([],[],linestyle='none', marker= item.get_marker(), markersize=item.get_markersize(), markeredgewidth=item.get_markeredgewidth(), markerfacecolor=item.get_markerfacecolor(), markeredgecolor=item.get_markeredgecolor())
    legend_item.set_linestyle('none')
    legend_items.append(legend_item)
legend_items.append(Line2D([],[],linestyle='none'))
legend_items.append(Line2D([],[],linestyle = '',  color='black', 
                           linewidth= .5,marker = 'd',
                           markerfacecolor='black',markeredgecolor = 'black',
                           markersize = 6, markeredgewidth = 1))
for linestyle in line_styles[:2]:
    legend_items.append(Line2D([],[],linestyle = linestyle, color='black', linewidth= 1,marker = 'd',
                               markerfacecolor='none', markeredgecolor = 'black',markersize = 2, markeredgewidth = 1))
entries = entries[:len(fig1B_data.columns)]
entries.append('')
entries.append('experiment')
for name, _ in main_paths:
    entries.append(name+' model')
# move legend outside of plot
ax.legend(legend_items,entries, **{'bbox_to_anchor':(1.05, 1), 'loc':2, 'borderaxespad':0.})

ax.set_xticks(ticks=[750,800,850,900,950])
ax.set_xlim(740,960)
ax.set_yticks(ticks=[-40,-30,-20,-10,0])
plt.savefig(os.path.join(image_path, '1b_by_conversion_sans_propene_all_models.pdf'), bbox_inches = 'tight')

################################
# figure 6 & table 2 - slopes of enrichment
print('creating figure 6 and table 2')
# simulate
model_to_slope_dict = {}
model_to_temperature_varying_enrichments = {}
for model_index in range(4):
    # get model data
    model_name = main_paths[model_index][0]
    model_path = main_paths[model_index][1]
    molecule_to_cluster_number = model_molecule_to_cluster_number[model_name]
    model_file = os.path.join(model_path, 'chem.cti')
    isotopomer_info = pd.read_csv(os.path.join(model_path,'isotopomer_cluster_info.csv'),index_col = 'name')

    # simulate model
    temperature_varying_enrichments = {}
    slopes_found = pd.DataFrame(index = ["dC2H4 = f(dCH4)",'dC2H6 = f(dCH4)','dC2H6 = f(dC2H4)','dBulk = f(dCH4)'])

    for temp in (800,850,900,950):
        delta_enrichments = pd.DataFrame(index=list(molecule_to_cluster_number.keys()) + ['bulk'])
        for psia in np.linspace(-10, 20, 5):
            solution = ctt.create_mechanism(model_file)

            center_labeled_delta = -28
            edge_labeled_delta = -28 + psia
            edge_labeled_fraction = am.getEnrichementFractionFromDelta(edge_labeled_delta)
            center_labeled_fraction = am.getEnrichementFractionFromDelta(center_labeled_delta)
            fraction_propane = 0.0049 # see supplemental
            initialMoleFractions_fig2={
                "CCC": fraction_propane * (1-center_labeled_fraction) * (1-edge_labeled_fraction)**2,
                "CCC-2": fraction_propane * center_labeled_fraction * edge_labeled_fraction**2,
                "CCC-3": fraction_propane * edge_labeled_fraction**2*(1-center_labeled_fraction),
                "CCC-4": fraction_propane * 2*edge_labeled_fraction * (1-edge_labeled_fraction) * center_labeled_fraction,
                "CCC-5": fraction_propane * 2*edge_labeled_fraction *(1-center_labeled_fraction) * (1-edge_labeled_fraction),
                "CCC-6": fraction_propane * center_labeled_fraction*(1-edge_labeled_fraction)**2,
                "[He]": 1-fraction_propane,
            }
            # set initial conditions of solution in kelvin pascals and mole fractions
            conditions = temp+273, 2e5, initialMoleFractions_fig2
            t_final = 95. / temp # see supplemental info

            output = ctt.run_simulation(solution, [0,t_final],conditions, 
                                  condition_type = 'constant-temperature-and-pressure',
                                  output_species = True,
                                  output_reactions = False,
                                  output_directional_reactions = False,
                                  output_rop_roc = False)
            species = output['species'].loc[t_final,:]
            isotopomer_info['conc'] = species

            # find enrichments
            enrichments_temp = {}
            for molecule, cluster_num in molecule_to_cluster_number.items():
                labels = isotopomer_info[isotopomer_info.cluster_number == cluster_num].index
                if molecule != 'H-atom':
                    enrichments_temp[molecule] = am.getDelta(isotopomer_info.loc[labels,'conc'],
                                               isotopomer_info.loc[labels,'enriched_atoms'],
                                               isotopomer_info.loc[labels,'unenriched_atoms'],
                                               )

            # add bulk term to values, which is taken at begining of simulation
            isotopomer_info['conc'] = output['species'].loc[0,:]
            labels = isotopomer_info[isotopomer_info.cluster_number == molecule_to_cluster_number['propane']].index
            enrichments_temp['bulk'] = am.getDelta(isotopomer_info.loc[labels,'conc'],
                                               isotopomer_info.loc[labels,'enriched_atoms'],
                                               isotopomer_info.loc[labels,'unenriched_atoms'],
                                               )
            # save data for next simulation
            delta_enrichments[psia] = pd.Series(enrichments_temp)

        #save delta_enrichements values
        delta_enrichments = delta_enrichments.T
        temperature_varying_enrichments[temp] = delta_enrichments
        # coordinates of the relationships found.
        axes_x = ['methane','methane','ethene','methane']
        axes_y = ['ethene','ethane','ethane','bulk']
        fits = []
        for index in range(len(axes_x)):
            try:
                fits.append(sm.ols(formula='{} ~ {}'.format(axes_y[index],axes_x[index]), data= delta_enrichments).fit())
            except:
                print(axes_x[index], axes_y[index])
                print(delta_enrichments)
                raise
        slopes = [] # list of slopes and intecept for each fit in fits
        for fit in fits:
            slope = fit.params.iat[1] # slope
            std_err_slope = fit.bse.iat[1] #uncertainty
            #if std_err_slope > abs(slope/1000):
            #    raise Exception('the standard error, {0}, is greater than one thousanth the slope, {1} for fit at temperature {2}.'.format(std_err_slope, slope, temp))
            slopes.append(slope)
        slopes_found[temp] = slopes
        # analyze for the slopes at each temperature.
    model_to_slope_dict[model_name] = slopes_found
    model_to_temperature_varying_enrichments[model_name] = temperature_varying_enrichments

# plot figure 6
delta_enrichments = model_to_temperature_varying_enrichments['full'][850]
axes_x = ['ethene','methane','methane','methane']
axes_y = ['ethane','ethane','ethene','bulk']
axes_x_ticks = [[-34,-28,-22,-16,-10],[-34,-28,-22,-16,-10],[-34,-28,-22,-16,-10], [-34,-28,-22,-16,-10]]
axes_y_ticks = [np.linspace(-48,-48+32,5),np.linspace(-48,-48+32,5),np.linspace(-40,-40+28,5),np.linspace(-40,-40+28,5)]
axes_y_lim = [(x.min(),x.max()) for x in axes_y_ticks]
axes_x_lim = [(min(x),max(x)) for x in axes_x_ticks]

sns.set_style('whitegrid')
sns.set_context('paper',font_scale=1.25)

f,ax = plt.subplots(ncols=2,nrows=2,sharex=False,sharey=False,)
f.subplots_adjust(hspace=.3,wspace=.5)
axes = ax.flatten()
index2section = {2:'A',1:'B',0:'C',3:'D'}
for index, ax in enumerate(axes):
    gilbert_data = pd.read_csv('../exp_data/Gilbert_fig2{}_engauge.csv'.format(index2section[index]))
    gilbert_data.plot(x=gilbert_data.columns[0], y=gilbert_data.columns[1], ax=ax,style = 'o')
    delta_enrichments.plot(x=axes_x[index], y=axes_y[index], ax=ax, style = '-')

    ax.set_ylabel("$\delta^{{13}}C_{{ {0} }}$\n$(\perthousand)$".format(axes_y[index]),rotation='horizontal',labelpad=25)
    ax.set_xlabel("$\delta^{{13}}C_{{ {0} }}(\perthousand)$".format(axes_x[index]))
    
    ax.set_xlim(axes_x_lim[index])
    ax.set_ylim(axes_y_lim[index])
    
    ax.set_xticks(axes_x_ticks[index])
    if index in [0,1]:
        ax.set_xticklabels(['','','','',''])
    ax.set_yticks(axes_y_ticks[index])
    if index in [1,3]:
        ax.set_yticklabels(['','','','',''])
    ax.legend([])

plt.savefig(os.path.join(image_path, '2.pdf'), bbox_inches = 'tight')

# get experimental data used in table 2
gilbert_data = pd.read_csv('../exp_data/Gilbert_table2_values.csv',index_col='relationship',)
gilbert_data.columns = [np.int64(temperature) for temperature in gilbert_data.columns]
gilbert_uncertainty = pd.read_csv('../exp_data/Gilbert_table2_uncertainty.csv',index_col='relationship',)

# get table describing scaled difference between model and experiment
scaled_differences = pd.DataFrame()
for model, slopes_found in model_to_slope_dict.items():
    difference = gilbert_data - slopes_found
    scaled_difference = difference.values.flatten() / gilbert_uncertainty.values.flatten()
    scaled_differences[model] = scaled_difference

scaled_differences['species_slopes'] = [value for value in gilbert_data.index.values for x in range(4)]
scaled_differences['temperature'] = [value for x in range(4) for value in gilbert_data.columns.values]
# get table of standard deviations between model and experiment
error_table = pd.DataFrame()
error_table['all'] = scaled_differences.std(axis='index')
for t in gilbert_data.columns:
    t1 = scaled_differences[scaled_differences.temperature == t]
    error_table[t] = t1.std(axis='index')
error_table = error_table.T
del error_table['temperature']
error_table = error_table[['full','drg','6rxn','3rxn']]
print('########### Table 2 ##############')
print(error_table.round(decimals=1).to_latex())

print('finished creating paper figures')
