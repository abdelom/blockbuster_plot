import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from matplotlib.ticker import FuncFormatter
import statsmodels.api as sm
import seaborn as sns
import demes
import demesdraw
import numpy as np
from scipy.stats import chi2


class Scenario:
    def __init__(self, likelihood, distance, times_ne, theta, se_thetas, residues, fitted_sfs):
        self.likelihood = likelihood  # La vraisemblance du SFS observé
        self.distance = distance  # La distance entre SFS observé et théorique
        self.times = times_ne  # Les indices associés aux temps
        self.theta = theta  # Les taux de mutation globaux

        # Nouveaux arguments
        self.se_thetas = se_thetas  # Les valeurs de theta (taux de mutation pour chaque population)
        self.residues = residues  # Les résidus entre le SFS observé et ajusté
        self.fitted_sfs = fitted_sfs  # Le SFS ajusté

        self.effective_sizes = []  # Les tailles de population effective (Ne)
        self.times_generations = []  # Les temps exprimés en générations
        self.times_years = []  # Les temps exprimés en années (si applicable)
        self.aic = self.calculate_aic()

    def calculate_aic(self):
        """
        Calcul de l'AIC (Akaike Information Criterion) pour ce scénario.
        L'AIC est donné par la formule : AIC = 2k - 2 * log(L)
        où k est le nombre de paramètres et L est la log-vraisemblance.
        """
        num_parameters = 2 * len(self.times)  # Breakpoints plus un paramètre
        aic_value = 2 * num_parameters - 2 * self.likelihood  # Formule de l'AIC
        return aic_value

    def __str__(self):
        """
        Retourne une chaîne lisible résumant le scénario.
        """
        text = [
            f"--- Scenario Summary ---",
            f"Likelihood       : {self.likelihood:.4f}",
            f"Distance         : {self.distance:.6f}",
            f"AIC              : {self.aic:.4f}",
            f"Theta (global)   : {self.theta}",
            f"Thetas (by pop)  : {self.se_thetas}",
            f"Times (Ne units) : {self.times}",
        ]
        if self.times_generations:
            text.append(f"Times (generations): {self.times_generations}")
        if self.times_years:
            text.append(f"Times (years): {self.times_years}")
        if self.effective_sizes:
            text.append(f"Effective sizes (Ne): {self.effective_sizes}")
        return "\n".join(text)



def parse_arguments():
    """Parse les arguments en ligne de commande et retourne le fichier d'entrée et le dossier de sortie."""
    parser = argparse.ArgumentParser(description="Parse a file and plot log likelihood vs number of breakpoints.")
    
    # Ajout des arguments d'entrée et de sortie
    parser.add_argument("-i", "--input", required=True, type=str, help="Path to the input file.")
    parser.add_argument("-o", "--output", required=True, type=str, help="Directory to save the plot.")
    
    # Nouveaux arguments optionnels
    parser.add_argument("-mu", "--mutation_rate", type=float, default=-1, help="Mutation rate per site per generation (default: -1)")
    parser.add_argument("-l", "--genome_length", type=float, default=-1, help="Genome length (default: -1)")
    parser.add_argument("-g", "--generation_time", type=float, default=-1, help="Generation time (default: -1)")
    
    # Nouvel argument pour les données de courbe constante par morceaux
    parser.add_argument("-pc", "--piecewise_data", nargs="+", type=float, default=None,
                        help="List of values for the piecewise constant curve (format: time1 time2 ... size1 size2 ...)")
    
    # Nouvel argument pour le facteur z pour la droite affine
    # parser.add_argument("-z", "--z_factor", type=float, default=None, help="Factor for affine line (default: None)")
    # parser.add_argument("-t", "--present_theta", type=float, default=None, help="theta in the present")
    parser.add_argument("-p", '--plot', action='store_false', help="outplot")
    parser.add_argument("-or", "--oriented", type=int, default=1, help="whether the sfs is oriented or not")
    args = parser.parse_args()
    
    # Vérifier si le fichier d'entrée existe
    if not os.path.exists(args.input):
        print("Error: Input file not found.")
        exit(1)
    
    # Vérifier si le dossier de sortie existe, sinon le créer
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    return args.input, args.output, args

def parse_document(file_path):
    # Liste pour stocker tous les scénarios
    scenarios = []
    n = None  # Variable pour stocker le nombre de points d'échelle de temps

    with open(file_path, "r") as file:
        # Lire chaque ligne du fichier
        for line in file:
            line = line.strip()
            
            # Si la ligne commence par ">", on commence à traiter un scénario
            if line.startswith(">"):
                # Initialiser les attributs du scénario
                likelihood = None
                distance = None
                breakpoints = []
                se_thetas = []
                residues = []
                fitted_sfs = []
                theta = None  # Valeur de theta globale (s'il est spécifié dans le fichier)
                
                # Parcourir les lignes suivantes pour lire les données associées au scénario
                while True:
                    next_line = file.readline().strip()
                    
                    if not next_line:  # Si ligne vide, fin du scénario
                        break
                    elif next_line.startswith("log_likelihood"):
                        likelihood = float(next_line.split(":")[-1].strip())
                    elif next_line.startswith("distance"):
                        distance = float(next_line.split(":")[-1].strip())
                    elif next_line.startswith("time"):
                        times_ne = list(map(float,next_line.strip().split()[6:]))
                    elif next_line.startswith("thetas"):
                        theta = list(map(float, next_line.strip().split()[1:]))
                    elif next_line.startswith("residues"):
                        residues = list(map(float, next_line.strip().split()[1:]))
                    elif next_line.startswith("fitted_sfs"):
                        fitted_sfs = list(map(float, next_line.strip().split()[1:]))
                        break
                    elif next_line.startswith("se_thetas"):  # Si theta est spécifié, l'extraire
                        se_thetas =  list(map(float, next_line.strip().split()[1:]))


                # Créer un objet Scenario avec les nouveaux arguments
                scenario = Scenario(
                    likelihood, distance, times_ne, theta, se_thetas, residues, fitted_sfs
                )
               
               # scenario.calculate_aic()  # Si applicable, calculer l'AI
                scenarios.append(scenario)

    return scenarios



def convert_timeandtheta(time_scale, scenario):
    new_times = [0.] 
    thetas = [scenario.theta[0]]
    b_index = 0
    relative_size = 1.
    timem1  = 0.
    
    # Loop over the time scale to calculate the new times and thetas
    for index, time in enumerate(time_scale):
        new_times.append(new_times[-1] + (time - timem1) * relative_size)
        thetas.append(scenario.theta[b_index])
        timem1 = time
        
        if b_index < len(scenario.breakpoints) and scenario.breakpoints[b_index] == index + 1:
            b_index += 1
            thetas.append(scenario.theta[b_index])
            new_times.append(new_times[-1])
            relative_size = scenario.theta[b_index] / scenario.theta[0]
    return new_times , thetas


def convert_scenario(scenario, mu, l, generation_time):
    if mu > 0 and l > 0:
        # Calculate effective population size at each breakpoint
        scenario.effective_sizes = [theta / (4 * mu * l) for theta in scenario.theta]
        
        # Times in generations scaled by Ne at present
        present_ne = scenario.effective_sizes[0]
        scenario.times_generations = [time * present_ne for time in scenario.times]
        
        # Times in years if generation_time is provided
        if generation_time > 0:
            scenario.times_years = [time * generation_time for time in scenario.times_generations]
    
    return scenario

def deme_format(scenario, output, mu=-1, g=-1, L=-1):

    # output = f"{len(scenario.times + 1)}_epochs.yml"
    times = scenario.times
    sizes = scenario.theta
    unit = "Ne generations"
    if mu > 0 and L > 0:
        times = scenario.times_generations
        sizes = scenario.effective_sizes
        unit = "generations"
        if g > 0:
            times = scenario.times_years
            unit = "years"
            generation_time = g  # à condition qu’il existe
    with open(output, 'w') as f:
        f.write(f"time_units: {unit}\n")
        if unit == "years" and g !=-1:
            f.write(f"generation_time: {generation_time}\n")
        else:
            f.write(f"generation_time: 1  \n")
        f.write(f"demes:\n")
        f.write(f"  - name: B\n")
        f.write(f"    start_time: .inf\n")
        f.write(f"    epochs:\n")
        i = len(sizes) - 1
        if len(times) > 0:
            f.write(f"      - start_size: {sizes[i]}\n")
            f.write(f"        end_time: {2.*times[-1]}\n")
        for time in times[::-1]: 
            f.write(f"      - start_size: {sizes[i]}\n")
            f.write(f"        end_time: {time}\n")
            i -= 1
        f.write(f"      - start_size: {sizes[0]}\n")
        f.write(f"        end_time: {0.}\n")
    return times, unit


def plot_demographic_scenarios4(scenarios, output_directory, best, o = 0, mu=-1, l=-1, g=-1):
    # Créer le sous-dossier "demes_format" s'il n'existe pas
    demes_dir = os.path.join(output_directory, "demes_format")
    os.makedirs(demes_dir, exist_ok=True)

    for index, scenario in enumerate(scenarios):
        # Nom de base pour les fichiers
        # print(f"{len(scenario.times) + 1}_epochs")
        base_name = f"{len(scenario.times) + 1}_epochs"

        # Chemin complet pour le fichier .yml
        yml_path = os.path.join(demes_dir, f"{base_name}.yml")

        # Génération du fichier .yml avec la fonction deme_format
        times, unit = deme_format(scenario, yml_path, mu, g, l)
        # Charger le graphe à partir du fichier YML
        graph = demes.load(yml_path)

        # Couleur rouge pour toutes les populations
        if index == best:
            col = "red"
        else:
            col = "blue"
        colours = {deme.name: col for deme in graph.demes}

        # Créer la figure avec 2 subplots (démographie + SFS)
        fig, axes = plt.subplots(2, 1, figsize=(8, 8))

        # --- Subplot 1 : Démographie ---
        ax1 = axes[0]
        demesdraw.size_history(graph, ax=ax1, inf_ratio=1e-4, log_size=True, colours=colours)
        ax1.set_xscale("log")
        if len(times) > 0:
            ax1.set_xlim(left=min(times)/5., right=2*max(times))
         # Ajuster le label Y selon l’unité
        if unit == "Ne generations":
            ylabel = r"$\theta = 4 \, \mu \, L \, N_e$"
        else:
            ylabel = r"$N_e$"
        ax1.set_ylabel(ylabel, rotation=90, labelpad=15, va='center')
        ax1.set_xlabel(f"Time in {unit}")
        ax1.set_title(f"Demographic history ({len(scenario.times) + 1} epochs)")

        # --- Subplot 2 : SFS observée et ajustée ---
        ax2 = axes[1]
        if hasattr(scenario, "fitted_sfs") and hasattr(scenario, "residues") \
           and scenario.fitted_sfs is not None and scenario.residues is not None:
            x = np.arange(1, len(scenario.fitted_sfs) + 1)
            fitted = np.array(scenario.fitted_sfs)
            observed = fitted + np.array(scenario.residues)

            # Multiplier chaque bin par son indice
            fitted_scaled = x * fitted
            observed_scaled = x * observed
            if o == 0:
                fitted_scaled[:-1] =  (2 * len(fitted_scaled) - x[:-1]) *  fitted_scaled[:-1] / (2 * len(fitted_scaled))
                observed_scaled[:-1] =  (2 * len(observed_scaled) - x[:-1]) *  observed_scaled[:-1] / (2 * len(observed_scaled))


     # Fitted SFS : courbe bleue continue
            ax2.plot(x, fitted_scaled, color=col, linewidth=1.5, label="Fitted SFS", zorder=3)
            # Observed SFS : petits points noirs transparents
            ax2.scatter(
                x, observed_scaled,
                color="black",
                s=20,              # petite taille
                alpha=0.6,         # transparence
                label="Observed SFS",
                zorder=2
            )
    
            # Mise en forme du subplot
            ax2.set_title("Observed vs Fitted Site Frequency Spectrum (SFS)")
            ax2.set_xlabel("Frequency class (i)")
            ax2.set_ylabel("i × E[Xi]")
            if o == 0:
                ax2.set_ylabel("i (n-i) × E[Xi] / n")

            ax2.legend()
            # ax2.grid(True, linestyle="--", alpha=0.4)
        else:
            ax2.text(0.5, 0.5, "No fitted SFS or residues available", ha="center", va="center")
            ax2.axis("off")

        # Ajustement des espacements
        plt.tight_layout()

        # Enregistrement de la figure
        fig_path = os.path.join(demes_dir, f"{base_name}.png")
        plt.savefig(fig_path, bbox_inches="tight", dpi=300)

        # Fermer la figure
        plt.close(fig)


def select_best_model_lrt_cumulative(scenarios, outputdirectory, alpha=0.05):
    """
    Sélectionne le meilleur modèle selon un test du rapport de vraisemblance cumulatif.
    
    Paramètres
    ----------
    scenarios : list
        Liste de modèles ordonnés par nombre d'epochs croissant.
        Chaque modèle doit avoir les attributs :
        - likelihood : log-vraisemblance
        - aic : AIC (optionnel)
    outputdirectory : str
        Chemin vers le dossier où sauvegarder le graphique.
    alpha : float
        Seuil de significativité pour le test du rapport de vraisemblance.
    show_plot : bool
        Si True, affiche la courbe des log-likelihoods.
    
    Retour
    ------
    best_index : int
        Indice du modèle retenu (Python index à partir de 0).
    """
    n_models = len(scenarios)
    epochs = np.arange(1, n_models + 1)
    logL = np.array([s.likelihood for s in scenarios])

    # --- Tracer et sauvegarder la courbe des log-likelihoods ---
    plt.figure(figsize=(6,4))
    plt.plot(epochs, logL, marker='o', color='blue', label='Log-likelihood')
    plt.xlabel("Number of epochs")
    plt.ylabel("Log-likelihood")
    plt.title("Model likelihoods vs number of epochs")
    plt.legend()
    plt.tight_layout()

    # Créer le dossier s'il n'existe pas
    os.makedirs(outputdirectory, exist_ok=True)
    plot_path = os.path.join(outputdirectory, "likelihood_vs_epochs.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()

    # --- Test du rapport de vraisemblance cumulatif ---
    best_index = 0  # modèle de référence : 1 epoch
    for i in range(1, n_models):
        ddl = 2 * (i - best_index)
        lr = 2 * (logL[i] - logL[best_index])
        critical = chi2.ppf(1 - alpha, df=ddl)
        if lr > critical:
            best_index = i

    print(f"✅ Best model according to cumulative LRT: {epochs[best_index]} epochs "
          f"(logL = {logL[best_index]:.3f})")
    print(f"📈 Plot saved in: {plot_path}")

    return best_index


# --- Fonctions pour chaque sous-plot ---

def plot_residuals_vs_fitted(ax, fitted, residuals):
    """Plot: Residuals vs Fitted Values"""
    ax.scatter(range(1, len(fitted) + 1), residuals, alpha=0.7)
    ax.axhline(0, color='red', linestyle='--', linewidth=1)
    ax.set_title('Residuals vs Fitted Values')
    ax.set_xlabel('Fitted Values')
    ax.set_ylabel('Residuals')


def plot_qq(ax, residuals):
    """Plot: Q-Q Plot"""
    sm.qqplot(residuals, line='45', fit=True, ax=ax)
    ax.set_title("Q-Q Plot of Residuals")


def plot_histogram(ax, residuals):
    """Plot: Residuals Histogram"""
    sns.histplot(residuals, kde=True, bins=20, ax=ax)
    ax.set_title('Residuals Histogram')
    ax.set_xlabel('Residuals')
    ax.set_ylabel('Density')


def plot_scale_location(ax, fitted, residuals):
    """Plot: Scale-Location Plot"""
    standardized_residuals = residuals / np.std(residuals)
    ax.scatter(fitted, np.sqrt(np.abs(standardized_residuals)), alpha=0.7)
    ax.axhline(0, color='red', linestyle='--', linewidth=1)
    ax.set_title('Scale-Location Plot')
    ax.set_xlabel('Fitted Values')
    ax.set_ylabel('√(|Standardized Residuals|)')


# --- Fonction globale pour le panneau ---

def plot_diagnostics(scenarios, output_directory):
    """
    Generates and saves diagnostic plots for a list of scenarios.
    
    Plots are saved in a subfolder 'diagnostic_plots' inside the output directory.
    Each plot is numbered according to the scenario order.
    
    Arguments:
    - scenarios: List of tuples (fitted, residuals), where fitted is a list of fitted values
                 and residuals is a list of residuals.
    - output_directory: Directory where the plots will be saved.
    """
    # Create a "diagnostic_plots" subfolder in the output directory
    diagnostics_dir = os.path.join(output_directory, "diagnostic_plots")
    os.makedirs(diagnostics_dir, exist_ok=True)

    # Compute the global y-axis limits for the top-left subplot (residuals vs fitted)
    all_residues = np.concatenate([np.array(s.residues) for s in scenarios])
    ymin, ymax = all_residues.min(), all_residues.max()

    # Loop over each scenario to create the figures
    for idx, scenario in enumerate(scenarios):
        # Create the figure with a 2x2 grid
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))
        fitted = np.array(scenario.fitted_sfs)
        residues = np.array(scenario.residues)

        # Plot the top-left subplot with shared y-axis
        plot_residuals_vs_fitted(axs[0, 0], fitted, residues)
        axs[0, 0].set_ylim(ymin, ymax)  # share the same y-axis across all scenarios

        # Plot the other subplots
        plot_qq(axs[0, 1], residues)
        plot_histogram(axs[1, 0], residues)
        plot_scale_location(axs[1, 1], fitted, residues)

        # Adjust spacing between subplots
        plt.tight_layout()

        # Save the figure in the subfolder
        output_file = os.path.join(diagnostics_dir, f"diagnostic_plot_{idx + 1}_epochs.png")
        plt.savefig(output_file)

        # Close the figure to free memory
        plt.close(fig)

    # Confirmation message
    print(f"All diagnostic plots have been saved in the folder: {diagnostics_dir}")


def main():
    # Parsing des arguments
    _, _, args = parse_arguments()
    
    # Parsing du fichier d'entrée
    scenarios = parse_document(args.input)   

    
    # Conversion des scénarios avec les paramètres mutation_rate, genome_length et generation_time
    for scenario in scenarios:
        convert_scenario(scenario, args.mutation_rate, args.genome_length, args.generation_time)
    
   # pdf_path = generate_pdf_report(scenarios, sfs, args.output)
   # print(f"PDF report saved at {pdf_path}")

    
    # Vérification et gestion de `piecewise_data` pour le tracé des courbes constantes par morceaux
    # import pdb; pdb.set_trace()
    # if args.piecewise_data is not None:
    #     piecewise_data = args.piecewise_data
    #     print("Piecewise data:", piecewise_data)
    # else:
    #     piecewise_data = None

    # # Vérification et gestion de `z_factor` pour la droite affine
    # if args.z_factor is not None:
    #     z_factor = args.z_factor
    #     print(f"Affine line factor z: {z_factor}")
    # else:
    #     z_factor = None
    # print(args.present_theta)
    # if args.plot:
    # Tracer les scénarios démographiques avec les courbes constantes par morceaux (si fournies)
    best = select_best_model_lrt_cumulative(scenarios, args.output, alpha=0.05)
    plot_demographic_scenarios4(
        scenarios, args.output, best, args.oriented,
        args.mutation_rate, args.genome_length, args.generation_time)
    plot_diagnostics(scenarios, args.output)
    print("Best model has", len(scenarios[best].times) + 1, "epochs")
        # # # Generate individual plots
        # plot_log_likelihood_vs_breakpoints(scenarios, os.path.join(args.output, "log_likelihood_plot.png"))
        # plot_aic_vs_breakpoints(scenarios, os.path.join(args.output, "aic_plot.png"))
    
    # Écriture des scénarios mis à jour dans un fichier de sortie
    # parsed_output_file = os.path.join(args.output, "parsed_output.txt")
    # write_scenario_output(
    #     parsed_output_file, time_scale, scenarios, 
    #     args.mutation_rate, args.genome_length, args.generation_time
    # )


if __name__ == "__main__":
    main()

def find_best_scenario(scenarios):
    """
    Finds the index of the scenario with the lowest AIC value.
    
    Parameters:
    - scenarios: List of Scenario objects
    
    Returns:
    - The index of the scenario with the lowest AIC value
    """
    if not scenarios:
        raise ValueError("The list of scenarios is empty.")
    scnearios2 = []
    for scenario in scenarios:
        flag = 1
        for theta in scenario.theta:
            if theta < 0:
                flag =  0
        if flag:
            scnearios2.append(scenario)
    # Find the index of the scenario with the lowest AIC
    best_index = min(range(len(scnearios2)), key=lambda i: scnearios2[i].aic)
    return best_index


# def find_best_scenariol(scenarios):
#     """
#     Finds the index of the scenario with the lowest AIC value.
    
#     Parameters:
#     - scenarios: List of Scenario objects
    
#     Returns:
#     - The index of the scenario with the lowest AIC value
#     """
#     if not scenarios:
#         raise ValueError("The list of scenarios is empty.")
#     for index, scenario in enumerate(scenarios):
#         if index == 0:
#             continue
#         if(2*(scenario.likelihood - scenarios[index - 1].likelihood) < 5.99):
#             return index - 1
#     # Find the index of the scenario with the lowest AIC
#     return len(scenarios) - 1

# def plot_log_likelihood_vs_breakpoints(scenarios, output_file):
#     """Trace la vraisemblance en fonction du nombre de breakpoints et enregistre le graphique dans un fichier."""
#     # Extraire les valeurs de vraisemblance et le nombre de breakpoints pour chaque scénario
#     log_likelihoods = [scenario.likelihood for scenario in scenarios]
#     num_breakpoints = [len(scenario.breakpoints) for scenario in scenarios]

#     # Tracer le graphe
#     plt.figure(figsize=(10, 6))
#     plt.plot(num_breakpoints[2:], log_likelihoods[2:], 'o-', color='blue')
#     plt.xlabel("Number of Breakpoints")
#     plt.ylabel("Log Likelihood")
#     plt.title("Log Likelihood vs Number of Breakpoints")
#     plt.grid(False)

#     # Sauvegarder le graphique dans le fichier spécifié
#     plt.savefig(output_file, dpi = 300)
#     print(f"Plot saved to {output_file}")


# def plot_aic_vs_breakpoints(scenarios, output_file):
#     """Trace l'AIC en fonction du nombre de breakpoints et enregistre le graphique dans un fichier."""
#     # Extraire les valeurs d'AIC et le nombre de breakpoints pour chaque scénario
#     aics = [scenario.aic for scenario in scenarios]  # Utilisation de l'attribut aic
#     num_breakpoints = [2 * len(scenario.breakpoints) for scenario in scenarios]  # Nombre de breakpoints pour chaque scénario
    
#     # Tracer le graphe
#     plt.figure(figsize=(10, 6))
#     plt.plot(num_breakpoints, aics, 'o-', color='red')
#     plt.xlabel("Number of Breakpoints")
#     plt.ylabel("AIC")
#     plt.title("AIC vs Number of parameters")
#     plt.grid(False)

#     # Sauvegarder le graphique dans le fichier spécifié
#     plt.savefig(output_file, dpi = 300)
#     print(f"Plot saved to {output_file}")

# def plot_demographic_scenarios3(scenarios, time_scale, output_directory, mu=-1, l=-1, g=-1, piecewise_data=None, z=None, theta=1e6):
#     """
#     Plot individual demographic scenarios, organize them on A4 pages with shared axis labels,
#     and save each individual plot as a separate PNG file. The first scenario is plotted in red.
#     """
#     # Create the main output directory and the "individual_plots" subdirectory
#     times1, times2 = [], []
#     os.makedirs(output_directory, exist_ok=True)
#     plots_directory = os.path.join(output_directory, "individual_plots")
#     os.makedirs(plots_directory, exist_ok=True)

#     # Prepare data for the combined plot if mu and l are provided
#     combined_data = []  # Will store (times, Nes, label) for each scenario
#     times = [[], [], []]
#     best = find_best_scenariol(scenarios)
#     # A4 page configuration
#     n_cols = 2  # Number of columns
#     n_rows = len(scenarios) // 2 + len(scenarios) % 2  # Number of rows per pag
#     if len(scenarios) > 1:
#         fig_a4, axs_a4 = plt.subplots(n_rows, n_cols, figsize=(8.27, 11.69), sharex=True, sharey=True)  # A4 size in inches
#         fig_a4.subplots_adjust(hspace=0.5, wspace=0.4)
#         axs_a4 = np.array(axs_a4)  # Ensure axs_a4 is a numpy array for easier indexing
#         plot_idx = 0  # Index for subplot tracking

#     for idx, scenario in enumerate(scenarios):
#         # Convert times and thetas for the scenario
#         times[0], thetas = convert_timeandtheta(time_scale, scenario)
#         flag = 0

#         for thetat in scenario.theta:
#             if thetat < 0:
#                 flag = 1
#         if flag:
#             continue
#         # Comp, time2 = [ute Ne if mu and l are provided
#         if mu != -1 and l != -1:
#             Nes = np.array(thetas) / (4 * mu * l)  # Effective population size
#             times[1] = np.array(times[0]) * Nes[0] *2 # Times in Ne generations
#             if g != -1:
#                 times[2] = times[1] * g  # Times in years
#         else:
#             Nes = None
#             times[1] = []  # No conversion if mu and l are not provided
#         # Add data to the combined plot if necessary
#         if Nes is not None:
#             combined_data.append((times[1] if g == -1 else times[2], Nes, f'Scenario {idx + 1}'))
#         else:
#             combined_data.append((times[0], thetas, f'Scenario {idx + 1}'))
#         scenario.times = garder_doublons(times[0])
#         if len(times[1]) !=0:
#             scenario.times_generations = garder_doublons(times[1])
#         if len(times[2]) != 0:
#             scenario.times_years = garder_doublons(times[2])

#         # # Plot the first scenario in red, and others in blue
#         color = 'red' if idx == best else 'blue'
#         fig_individual, ax_individual = plt.subplots(figsize=(6, 4))  # Individual plot size
#         ax_individual.set_xscale('log')
#         ax_individual.set_yscale('log')
#         plot_individual_scenario(ax_individual, times, thetas, idx, Nes, piecewise_data, z, color, theta)
#         ax_individual.set_title(f"{idx + 1} epochs")
#         individual_output_file = os.path.join(plots_directory, f'scenario_{idx + 1}_epochs.png')
#         plt.tight_layout()
#         fig_individual.savefig(individual_output_file, dpi=300)
#         plt.close(fig_individual)

#     deme_format(scenarios[best], output_directory + "/deme.yml", g)
#     graph = demes.load(output_directory + "/deme.yml")
#     colours = {deme.name: "red" for deme in graph.demes}

#     # Tracer le graphe avec les bonnes couleurs
#     demesdraw.size_history(graph, log_time=True, colours=colours)
#     plt.show()
#     print(f"All individual plots have been saved in the directory: {plots_directory}")

# def write_scenario_output(output_file, time_scale, scenarios, mu, l, generation_time):
#     with open(output_file, 'w') as f:
#         # Write the original time scale
#         # Write each scenario with its converted expressions
#         for i, scenario in enumerate(scenarios):
#             # Write likelihood and distance
#             time, _ = convert_timeandtheta(time_scale, scenario)
#             scenario.times =garder_doublons(time)
#             f.write(f"> {i + 1} epochs\n")
#             f.write(f"lik : {scenario.likelihood} dist : {scenario.distance} aic : {scenario.aic}\n")
            
#             # Line 1: Breakpoints in coalescent units (original input format)
#             #f.write(" ".join(map(str, scenario.breakpoints)) + " ")
#             f.write("Populational mutation rate (theta): " + " ".join(f"{theta:.6f}" for theta in scenario.theta) + " ")
#             f.write("Times of change in Ne generations: " + " ".join(f"{time:.6f}" for time in scenario.times) + " \n")
#             if mu > 0 and l > 0:
#                 # Line 2: Effective population size (Ne) and times in generations
#                 #f.write(" ".join(map(str, scenario.breakpoints)) + " ")
#                 f.write("Effective size (Ne): " + " ".join(f"{ne:.6f}" for ne in scenario.effective_sizes) + " ")
#                 f.write("Times of change in generations: " + " ".join(f"{time_gen:.6f}" for time_gen in scenario.times_generations) + "\n")
                
#                 # Line 3: Effective population size (Ne) and times in years (if generation time is provided)
#                 if generation_time > 0:
#                     #f.write(" ".join(map(str, scenario.breakpoints)) + " ")
#                     f.write("Effective size (Ne): " + " ".join(f"{ne:.6f}" for ne in scenario.effective_sizes) + " ")
#                     f.write("Times of change in years: " + " ".join(f"{time_year:.6f}" for time_year in scenario.times_years) + "\n")

# def garder_doublons(liste):
#     result = []
#     i = 0
#     while i < len(liste):
#         # Vérifier si la valeur actuelle a un doublon consécutif
#         if i + 1 < len(liste) and liste[i] == liste[i + 1]:
#             result.append(liste[i])
#             # Passer à la fin des doublons consécutifs
#             while i + 1 < len(liste) and liste[i] == liste[i + 1]:
#                 i += 1
#         i += 1
#     return result

# def plot_individual_scenario(ax1, times, thetas, scenario_idx, Nes=None, piecewise_data=None, z=None, color="blue", theta=2e7):
#     """
#     Plot the individual scenario: either Ne or Theta values with time.
#     Adjust axis label sizes based on whether the plot is for a PDF or individual PNG.

#     Parameters:
#     - ax1: Matplotlib axis to plot on.
#     - times: List of time scales (Ne generations, generations, years).
#     - thetas: List of theta values.
#     - scenario_idx: Index of the scenario being plotted.
#     - Nes: Effective population sizes (if available).
#     - piecewise_data: Data for piecewise constant curves (optional).
#     - z: Parameters for affine function (optional).
#     - color: Color of the plot.
#     - is_pdf: Boolean, True if the plot is for a PDF, False otherwise.
#     """
#     # Adjust font sizes for PDF or individual plots
#     label_fontsize = 12
#     tick_fontsize = 10 
#     if Nes is not None:  # If Ne is calculated, plot it
#         if len(times[2]) > 0:  # If time is in years
#             ax1.plot(times[2], Nes, label=f'Scenario {scenario_idx + 1} (years)', color=color)
#             ax1.set_xlabel('Time (years)', fontsize=label_fontsize, fontweight='bold')
#         else:  # If time is in Ne generations
#             ax1.plot(times[1], Nes, color=color)
#             ax1.set_xlabel('Time (number of generations)', fontsize=label_fontsize, fontweight='bold')
#         ax1.set_ylabel('Effective Population \n Size (Ne)', fontsize=label_fontsize, fontweight='bold')
#     else:  # If Ne is not calculated, plot Theta
#         ax1.plot(times[0], thetas, label=f'Scenario {scenario_idx + 1}', color=color)
#         ax1.set_xlabel('Time (Ne generations)', fontsize=label_fontsize, fontweight='bold')
#         ax1.set_ylabel('Population Mutation \n Rate (Theta)', fontsize=label_fontsize, fontweight='bold')

#     # Handle piecewise constant curve if data is provided
#     if piecewise_data is not None:
#         plot_piecewise(ax1, piecewise_data, theta)

#     # Handle affine function plot if z is provided
#     if z is not None:
#         plot_affine_function(ax1, times, theta, z)

#     # Customize plot appearance
#     #     y_min = min(thetas+sizes_scaled) * 0.9  # Add some margin
#     # y_max = max(thetas+ sizes_scaled) * 1.1  # Add some margin
#     # ax1.set_ylim(y_min, y_max)
#     ax1.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:g}'))
#     ax1.tick_params(axis='x', labelsize=tick_fontsize, rotation=45)
#     ax1.tick_params(axis='y', labelsize=tick_fontsize)
#     ax1.grid(True, which='both', linestyle='', linewidth=0.5, alpha=0.7)
#     ax1.spines['top'].set_visible(False)
#     ax1.spines['right'].set_visible(False)


# def plot_piecewise(ax, piecewise_data, theta):
#     """
#     Plot a piecewise-constant population size curve using ax.plot
#     and remove top/right spines.
#     """
#     n = len(piecewise_data) // 2
#     times_piecewise = [0] + piecewise_data[:n] + [piecewise_data[:n][-1] + 1]
#     sizes_piecewise = [1] + piecewise_data[n:]
#     sizes_scaled = [s * theta for s in sizes_piecewise]

#     for i in range(len(sizes_scaled)):
#         # Horizontal segment
#         x_horiz = [times_piecewise[i], times_piecewise[i + 1]]
#         y_horiz = [sizes_scaled[i], sizes_scaled[i]]
#         ax.plot(x_horiz, y_horiz, color='black', label="Inferred" if i == 0 else "")

#         if i < len(sizes_scaled) - 1:
#             # Vertical segment
#             x_vert = [times_piecewise[i + 1], times_piecewise[i + 1]]
#             y_vert = [sizes_scaled[i], sizes_scaled[i + 1]]
#             ax.plot(x_vert, y_vert, color='black')

# def plot_affine_function(ax1, times, theta, z):
#     """
#     Plot the affine function if z is provided.
#     """
#     affine_y = theta * z * np.array(times[0]) + theta
#     ax1.plot(times[0], affine_y, label=f'Simulation', color='black')

