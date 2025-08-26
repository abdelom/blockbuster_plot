import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from matplotlib.ticker import FuncFormatter
import statsmodels.api as sm
import seaborn as sns


class Scenario:
    def __init__(self, likelihood, distance, breakpoints, theta, se_thetas, residues, fitted_sfs):
        self.likelihood = likelihood  # La vraisemblance du SFS observé
        self.distance = distance  # La distance entre SFS observé et théorique
        self.breakpoints = breakpoints  # Les indices associés aux temps
        self.theta = theta  # Les taux de mutation globaux

        # Nouveaux arguments
        self.se_thetas = se_thetas  # Les valeurs de theta (taux de mutation pour chaque population)
        self.residues = residues  # Les résidus entre le SFS observé et ajusté
        self.fitted_sfs = fitted_sfs  # Le SFS ajusté

        # Attributs existants
        self.times = []  # Les temps de changement de taille de population
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
        num_parameters = 2 * len(self.breakpoints)  # Breakpoints plus un paramètre
        aic_value = 2 * num_parameters - 2 * self.likelihood  # Formule de l'AIC
        return aic_value



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
    parser.add_argument("-z", "--z_factor", type=float, default=None, help="Factor for affine line (default: None)")
    parser.add_argument("-t", "--present_theta", type=float, default=None, help="theta in the present")
    parser.add_argument("-p", '--plot', action='store_false', help="outplot")
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
    # Liste pour stocker les valeurs de l'échelle de temps
    time_scale = []
    # Liste pour stocker tous les scénarios
    scenarios = []
    n = None  # Variable pour stocker le nombre de points d'échelle de temps

    with open(file_path, "r") as file:
        # Lire chaque ligne du fichier
        for line in file:
            line = line.strip()
            # Si la ligne n'est pas vide et ne commence pas par ">", elle fait partie de l'échelle de temps
            if line and not line.startswith(">"):
                time_scale.append(float(line))
            
            # Si la ligne commence par ">", on commence à traiter un scénario
            elif line.startswith(">"):
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
                    elif next_line.startswith("breakpoints"):
                        breakpoints = list(map(float,next_line.strip().split()[1:]))
                    elif next_line.startswith("thetas"):
                        theta = list(map(float, next_line.strip().split()[1:]))
                    elif next_line.startswith("residues"):
                        residues = list(map(float, next_line.strip().split()[1:]))
                    elif next_line.startswith("fitted_sfs"):
                        fitted_sfs = list(map(float, next_line.strip().split()[1:]))
                        break
                    elif next_line.startswith("se_thetas"):  # Si theta est spécifié, l'extraire
                        se_thetas =  list(map(float, next_line.strip().split()[1:]))

                # Calculer n si nécessaire
                if n is None:
                    n = len(time_scale)

                # Créer un objet Scenario avec les nouveaux arguments
                scenario = Scenario(
                    likelihood, distance, breakpoints, theta, se_thetas, residues, fitted_sfs
                )
               
               # scenario.calculate_aic()  # Si applicable, calculer l'AI
                scenarios.append(scenario)

    return time_scale, scenarios

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


def find_best_scenariol(scenarios):
    """
    Finds the index of the scenario with the lowest AIC value.
    
    Parameters:
    - scenarios: List of Scenario objects
    
    Returns:
    - The index of the scenario with the lowest AIC value
    """
    if not scenarios:
        raise ValueError("The list of scenarios is empty.")
    for index, scenario in enumerate(scenarios):
        if index == 0:
            continue
        if(2*(scenario.likelihood - scenarios[index - 1].likelihood) < 5.99):
            return index - 1
    # Find the index of the scenario with the lowest AIC
    return len(scenarios) - 1

def plot_log_likelihood_vs_breakpoints(scenarios, output_file):
    """Trace la vraisemblance en fonction du nombre de breakpoints et enregistre le graphique dans un fichier."""
    # Extraire les valeurs de vraisemblance et le nombre de breakpoints pour chaque scénario
    log_likelihoods = [scenario.likelihood for scenario in scenarios]
    num_breakpoints = [len(scenario.breakpoints) for scenario in scenarios]

    # Tracer le graphe
    plt.figure(figsize=(10, 6))
    plt.plot(num_breakpoints, log_likelihoods, 'o-', color='blue')
    plt.xlabel("Number of Breakpoints")
    plt.ylabel("Log Likelihood")
    plt.title("Log Likelihood vs Number of Breakpoints")
    plt.grid(True)

    # Sauvegarder le graphique dans le fichier spécifié
    plt.savefig(output_file, dpi = 300)
    print(f"Plot saved to {output_file}")


def plot_aic_vs_breakpoints(scenarios, output_file):
    """Trace l'AIC en fonction du nombre de breakpoints et enregistre le graphique dans un fichier."""
    # Extraire les valeurs d'AIC et le nombre de breakpoints pour chaque scénario
    aics = [scenario.aic for scenario in scenarios]  # Utilisation de l'attribut aic
    num_breakpoints = [2 * len(scenario.breakpoints) for scenario in scenarios]  # Nombre de breakpoints pour chaque scénario
    
    # Tracer le graphe
    plt.figure(figsize=(10, 6))
    plt.plot(num_breakpoints, aics, 'o-', color='red')
    plt.xlabel("Number of Breakpoints")
    plt.ylabel("AIC")
    plt.title("AIC vs Number of parameters")
    plt.grid(True)

    # Sauvegarder le graphique dans le fichier spécifié
    plt.savefig(output_file, dpi = 300)
    print(f"Plot saved to {output_file}")



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
        scenario.times_generations = [20 * time * present_ne for time in scenario.times]
        
        # Times in years if generation_time is provided
        if generation_time > 0:
            scenario.times_years = [time * generation_time for time in scenario.times_generations]
    
    return scenario

def garder_doublons(liste):
    result = []
    i = 0
    while i < len(liste):
        # Vérifier si la valeur actuelle a un doublon consécutif
        if i + 1 < len(liste) and liste[i] == liste[i + 1]:
            result.append(liste[i])
            # Passer à la fin des doublons consécutifs
            while i + 1 < len(liste) and liste[i] == liste[i + 1]:
                i += 1
        i += 1
    return result

def plot_individual_scenario(ax1, times, thetas, scenario_idx, Nes=None, piecewise_data=None, z=None, color="blue", theta=2e7):
    """
    Plot the individual scenario: either Ne or Theta values with time.
    Adjust axis label sizes based on whether the plot is for a PDF or individual PNG.

    Parameters:
    - ax1: Matplotlib axis to plot on.
    - times: List of time scales (Ne generations, generations, years).
    - thetas: List of theta values.
    - scenario_idx: Index of the scenario being plotted.
    - Nes: Effective population sizes (if available).
    - piecewise_data: Data for piecewise constant curves (optional).
    - z: Parameters for affine function (optional).
    - color: Color of the plot.
    - is_pdf: Boolean, True if the plot is for a PDF, False otherwise.
    """
    # Adjust font sizes for PDF or individual plots
    label_fontsize = 12
    tick_fontsize = 10 
    if Nes is not None:  # If Ne is calculated, plot it
        if len(times[2]) > 0:  # If time is in years
            ax1.plot(times[2], Nes, label=f'Scenario {scenario_idx + 1} (years)', color=color)
            ax1.set_xlabel('Time (years)', fontsize=label_fontsize, fontweight='bold')
        else:  # If time is in Ne generations
            ax1.plot(times[1], Nes, color=color)
            ax1.set_xlabel('Time (number of generations)', fontsize=label_fontsize, fontweight='bold')
        ax1.set_ylabel('Effective Population \n Size (Ne)', fontsize=label_fontsize, fontweight='bold')
    else:  # If Ne is not calculated, plot Theta
        ax1.plot(times[0], thetas, label=f'Scenario {scenario_idx + 1}', color=color)
        ax1.set_xlabel('Time (Ne generations)', fontsize=label_fontsize, fontweight='bold')
        ax1.set_ylabel('Population Mutation \n Rate (Theta)', fontsize=label_fontsize, fontweight='bold')

    # Handle piecewise constant curve if data is provided
    if piecewise_data is not None:
        plot_piecewise(ax1, piecewise_data, theta)

    # Handle affine function plot if z is provided
    if z is not None:
        plot_affine_function(ax1, times, theta, z)

    # Customize plot appearance
    #     y_min = min(thetas+sizes_scaled) * 0.9  # Add some margin
    # y_max = max(thetas+ sizes_scaled) * 1.1  # Add some margin
    # ax1.set_ylim(y_min, y_max)
    ax1.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:g}'))
    ax1.tick_params(axis='x', labelsize=tick_fontsize, rotation=45)
    ax1.tick_params(axis='y', labelsize=tick_fontsize)
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

def plot_piecewise(ax1, piecewise_data, theta):
    """
    Plot the piecewise constant curve (horizontal and vertical lines)
    using only ax.plot instead of hlines and vlines.
    """
    n = len(piecewise_data) // 2
    times_piecewise = [0] + piecewise_data[:n] + [piecewise_data[:n][-1] + 1]
    sizes_piecewise = [1] + piecewise_data[n:]

    # Scale the sizes
    sizes_scaled = [s * theta for s in sizes_piecewise]

    for i in range(len(sizes_scaled)):
        # Horizontal segment
        x_horiz = [times_piecewise[i], times_piecewise[i + 1]]
        y_horiz = [sizes_scaled[i], sizes_scaled[i]]
        ax1.plot(x_horiz, y_horiz, color='black', label="Inferred" if i == 0 else "")

        if i < len(sizes_scaled) - 1:
            # Vertical segment
            x_vert = [times_piecewise[i + 1], times_piecewise[i + 1]]
            y_vert = [sizes_scaled[i], sizes_scaled[i + 1]]
            ax1.plot(x_vert, y_vert, color='black')


def plot_affine_function(ax1, times, theta, z):
    """
    Plot the affine function if z is provided.
    """
    affine_y = theta * z * np.array(times[0]) + theta
    ax1.plot(times[0], affine_y, label=f'Simulation', color='black')


def deme_format(scenario, output):
    times = scenario.times
    sizes = scenario.theta
    unit = "Ne generations"
    if len(scenario.times_generations) > 0:
        times = scenario.times_generations
        sizes = scenario.effective_sizes
        unit = "generations"
    if len(scenario.times_years) > 0:
        times = scenario.times_generations
        unit = "years"
    with open(output, 'w') as f:
        f.write(f"time_units: {unit}\n")
        f.write(f"demes:\n")
        f.write(f"  - name: B\n")
        f.write(f"    start_time: .inf\n")
        f.write(f"    epochs:\n")
        i = len(sizes) - 1
        for times in times[::-1]: 
            f.write(f"      - start_size: {sizes[i]}\n")
            f.write(f"        end_time: {times}\n")
            i -= 1
        f.write(f"      - start_size: {sizes[0]}\n")
        f.write(f"        end_time: {0}\n")

def plot_demographic_scenarios3(scenarios, time_scale, output_directory, mu=-1, l=-1, g=-1, piecewise_data=None, z=None, theta=1e6):
    """
    Plot individual demographic scenarios, organize them on A4 pages with shared axis labels,
    and save each individual plot as a separate PNG file. The first scenario is plotted in red.
    """
    # Create the main output directory and the "individual_plots" subdirectory
    times1, times2 = [], []
    os.makedirs(output_directory, exist_ok=True)
    plots_directory = os.path.join(output_directory, "individual_plots")
    os.makedirs(plots_directory, exist_ok=True)

    # Prepare data for the combined plot if mu and l are provided
    combined_data = []  # Will store (times, Nes, label) for each scenario
    times = [[], [], []]
    best = find_best_scenariol(scenarios)
    # A4 page configuration
    n_cols = 2  # Number of columns
    n_rows = len(scenarios) // 2 + len(scenarios) % 2  # Number of rows per pag
    if len(scenarios) > 1:
        fig_a4, axs_a4 = plt.subplots(n_rows, n_cols, figsize=(8.27, 11.69), sharex=True, sharey=True)  # A4 size in inches
        fig_a4.subplots_adjust(hspace=0.5, wspace=0.4)
        axs_a4 = np.array(axs_a4)  # Ensure axs_a4 is a numpy array for easier indexing
        plot_idx = 0  # Index for subplot tracking

    for idx, scenario in enumerate(scenarios):
        # Convert times and thetas for the scenario
        times[0], thetas = convert_timeandtheta(time_scale, scenario)
        flag = 0

        for thetat in scenario.theta:
            if thetat < 0:
                flag = 1
        if flag:
            continue
        # Comp, time2 = [ute Ne if mu and l are provided
        if mu != -1 and l != -1:
            Nes = np.array(thetas) / (4 * mu * l)  # Effective population size
            times[1] = np.array(times[0]) * Nes[0] *2 # Times in Ne generations
            if g != -1:
                times[2] = times[1] * g  # Times in years
        else:
            Nes = None
            times[1] = []  # No conversion if mu and l are not provided
        # Add data to the combined plot if necessary
        if Nes is not None:
            combined_data.append((times[1] if g == -1 else times[2], Nes, f'Scenario {idx + 1}'))
        else:
            combined_data.append((times[0], thetas, f'Scenario {idx + 1}'))
        scenario.times = garder_doublons(times[0])
        if len(times[1]) !=0:
            scenario.times_generations = garder_doublons(times[1])
        if len(times[2]) != 0:
            scenario.times_years = garder_doublons(times[2])

        # Get current subplot for A4 layout
        # row, col = divmod(plot_idx, n_cols)
        # ax = axs_a4[row, col]

        # # Plot individual scenario in the subplot
        # ax.set_xscale('log')
        # ax.set_yscale('log')

        # # Plot the first scenario in red, and others in blue
        color = 'red' if idx == best else 'blue'

        # # Appel à plot_individual_scenario avec les bons arguments
        # plot_individual_scenario(ax, times, thetas, idx, Nes, piecewise_data, z, color, theta)

        # ax.set_title(f"{idx + 1} bloks")

        # # Show labels only for the first column (y-axis) and last row (x-axis)
        # ax.set_ylabel("")  # Remove y-axis label for non-first columns
        
        # ax.set_xlabel("")  # Remove x-axis label for non-last rows
        # ax.set_xlabel("")
        # ax.tick_params(axis='x', labelsize=15, rotation = 45)
        # ax.tick_params(axis='y', labelsize=15)
        # Save the individual plot as PNG
        fig_individual, ax_individual = plt.subplots(figsize=(6, 4))  # Individual plot size
        ax_individual.set_xscale('log')
        ax_individual.set_yscale('log')
        plot_individual_scenario(ax_individual, times, thetas, idx, Nes, piecewise_data, z, color, theta)
        ax_individual.set_title(f"{idx + 1} blocks")
        # ax_individual.set_xlabel("Time (log scale)", fontsize=15)
        # ax_individual.set_ylabel("Population Size (log scale)", fontsize=15)
        individual_output_file = os.path.join(plots_directory, f'scenario_{idx + 1}.png')
        plt.tight_layout()
        fig_individual.savefig(individual_output_file, dpi=300)
        plt.close(fig_individual)

        # plot_idx += 1

    #     # If the current page is full, save it and start a new one
    #     if plot_idx >= n_rows * n_cols:
    #         # Add shared labels before saving
    #         fig_a4.supxlabel("Time (log scale)", fontsize=6)  # Reduced font size for x-axis
    #         fig_a4.supylabel("Population Size (log scale)", fontsize=6)  # Reduced font size for y-axis
    #         a4_output_file = os.path.join(plots_directory, f'individual_plots_page_{idx // (n_rows * n_cols) + 1}.pdf')
    #         plt.tight_layout()
    #         fig_a4.savefig(a4_output_file)
    #         plt.close(fig_a4)

    #         # Create a new A4 figure for the next set of plots
    #         fig_a4, axs_a4 = plt.subplots(n_rows, n_cols, figsize=(8.27, 11.69), sharex=True, sharey=True)
    #         fig_a4.subplots_adjust(hspace=0.5, wspace=0.4)
    #         axs_a4 = np.array(axs_a4)  # Ensure axs_a4 is a numpy array
    #         plot_idx = 0

    # # Save the last page if it has any plots
    # if plot_idx > 0:
    #     # Add shared labels before saving
    #     fig_a4.supxlabel("Time (log scale)", fontsize=6)  # Reduced font size for x-axis
    #     fig_a4.supylabel("Population Size (log scale)", fontsize=6)  # Reduced font size for y-axis
    #     a4_output_file = os.path.join(plots_directory, f'individual_plots_page_{len(scenarios) // (n_rows * n_cols) + 1}.pdf')
    #     plt.tight_layout()
    #     fig_a4.savefig(a4_output_file)
    #     plt.close(fig_a4)
    deme_format(scenarios[best], output_directory + "/deme.yml")
    print(f"All individual plots have been saved in the directory: {plots_directory}")



# --- Fonctions pour chaque sous-plot ---

def plot_residuals_vs_fitted(ax, fitted, residuals):
    """Graphique : Résidus vs Valeurs ajustées"""
    ax.scatter(range(1, len(fitted) + 1), residuals, alpha=0.7)
    ax.axhline(0, color='red', linestyle='--', linewidth=1)
    ax.set_title('Résidus vs Valeurs ajustées')
    ax.set_xlabel('Valeurs ajustées')
    ax.set_ylabel('Résidus')

def plot_qq(ax, residuals):
    """Graphique : Q-Q Plot"""
    sm.qqplot(residuals, line='45', fit=True, ax=ax)
    ax.set_title("Q-Q plot des résidus")

def plot_histogram(ax, residuals):
    """Graphique : Histogramme des résidus"""
    sns.histplot(residuals, kde=True, bins=20, ax=ax)
    ax.set_title('Histogramme des résidus')
    ax.set_xlabel('Résidus')
    ax.set_ylabel('Densité')

def plot_scale_location(ax, fitted, residuals):
    """Graphique : Scale-Location Plot"""
    standardized_residuals = residuals / np.std(residuals)
    ax.scatter(fitted, np.sqrt(np.abs(standardized_residuals)), alpha=0.7)
    ax.axhline(0, color='red', linestyle='--', linewidth=1)
    ax.set_title('Scale-Location Plot')
    ax.set_xlabel('Valeurs ajustées')
    ax.set_ylabel('√(|Résidus standardisés|)')

# --- Fonction globale pour le panneau ---

def plot_diagnostics(scenarios, output_directory):
    """
    Génère et sauvegarde des graphiques de diagnostic pour une liste de scénarios.
    
    Les graphiques sont sauvegardés dans un sous-dossier 'diagnostic_plots' du répertoire de sortie.
    Chaque graphique est numéroté en fonction de l'ordre des scénarios.
    
    Arguments :
    - scenarios : Liste de tuples (fitted, residuals), où fitted est une liste des valeurs ajustées 
                  et residuals est une liste des résidus.
    - output_directory : Répertoire dans lequel les graphiques seront sauvegardés.
    """
    # Créer un sous-dossier "diagnostic_plots" dans le répertoire de sortie
    diagnostics_dir = os.path.join(output_directory, "diagnostic_plots")
    os.makedirs(diagnostics_dir, exist_ok=True)

    for idx, scenario in enumerate(scenarios):
        # Créer la figure pour les diagnostics
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))  # 2x2 grilles
        fitted = np.array(scenario.fitted_sfs)
        residues = np.array(scenario.residues)
        # Appeler chaque fonction pour tracer les sous-graphiques
        plot_residuals_vs_fitted(axs[0, 0], fitted, residues)
        plot_qq(axs[0, 1], residues)
        plot_histogram(axs[1, 0], residues)
        plot_scale_location(axs[1, 1], fitted, residues)

        # Ajuster l'espacement entre les graphiques
        plt.tight_layout()

        # Sauvegarder la figure dans le sous-dossier
        output_file = os.path.join(diagnostics_dir, f"diagnostic_plot_{idx + 1}.png")
        plt.savefig(output_file)

        # Fermer la figure pour libérer la mémoire
        plt.close(fig)

    print(f"Tous les diagnostics ont été sauvegardés dans le dossier : {diagnostics_dir}")


def write_scenario_output(output_file, time_scale, scenarios, mu, l, generation_time):
    with open(output_file, 'w') as f:
        # Write the original time scale
        # Write each scenario with its converted expressions
        for i, scenario in enumerate(scenarios):
            # Write likelihood and distance
            time, _ = convert_timeandtheta(time_scale, scenario)
            scenario.times =garder_doublons(time)
            f.write(f"> {i + 1} blocks\n")
            f.write(f"lik : {scenario.likelihood} dist : {scenario.distance} aic : {scenario.aic}\n")
            
            # Line 1: Breakpoints in coalescent units (original input format)
            #f.write(" ".join(map(str, scenario.breakpoints)) + " ")
            f.write("Populational mutation rate (theta): " + " ".join(f"{theta:.6f}" for theta in scenario.theta) + " ")
            f.write("Times of change in Ne generations: " + " ".join(f"{time:.6f}" for time in scenario.times) + " \n")
            if mu > 0 and l > 0:
                # Line 2: Effective population size (Ne) and times in generations
                #f.write(" ".join(map(str, scenario.breakpoints)) + " ")
                f.write("Effective size (Ne): " + " ".join(f"{ne:.6f}" for ne in scenario.effective_sizes) + " ")
                f.write("Times of change in generations: " + " ".join(f"{time_gen:.6f}" for time_gen in scenario.times_generations) + "\n")
                
                # Line 3: Effective population size (Ne) and times in years (if generation time is provided)
                if generation_time > 0:
                    #f.write(" ".join(map(str, scenario.breakpoints)) + " ")
                    f.write("Effective size (Ne): " + " ".join(f"{ne:.6f}" for ne in scenario.effective_sizes) + " ")
                    f.write("Times of change in years: " + " ".join(f"{time_year:.6f}" for time_year in scenario.times_years) + "\n")



def main():
    # Parsing des arguments
    _, _, args = parse_arguments()
    
    # Parsing du fichier d'entrée
    time_scale, scenarios = parse_document(args.input)
    
    # Conversion des scénarios avec les paramètres mutation_rate, genome_length et generation_time
    for scenario in scenarios:
        convert_scenario(scenario, args.mutation_rate, args.genome_length, args.generation_time)
    
   # pdf_path = generate_pdf_report(scenarios, sfs, args.output)
   # print(f"PDF report saved at {pdf_path}")

    
    # Vérification et gestion de `piecewise_data` pour le tracé des courbes constantes par morceaux
    if args.piecewise_data is not None:
        piecewise_data = args.piecewise_data
        print("Piecewise data:", piecewise_data)
    else:
        piecewise_data = None

    # Vérification et gestion de `z_factor` pour la droite affine
    if args.z_factor is not None:
        z_factor = args.z_factor
        print(f"Affine line factor z: {z_factor}")
    else:
        z_factor = None
    print(args.present_theta)
    if args.plot:
    # Tracer les scénarios démographiques avec les courbes constantes par morceaux (si fournies)
        plot_demographic_scenarios3(
            scenarios, time_scale, args.output, 
            args.mutation_rate, args.genome_length, args.generation_time, 
            piecewise_data=piecewise_data, z=z_factor, theta = args.present_theta
        )
        plot_diagnostics(scenarios, args.output)

        # # Generate individual plots
        plot_log_likelihood_vs_breakpoints(scenarios, os.path.join(args.output, "log_likelihood_plot.png"))
        plot_aic_vs_breakpoints(scenarios, os.path.join(args.output, "aic_plot.png"))
    
    # Écriture des scénarios mis à jour dans un fichier de sortie
    parsed_output_file = os.path.join(args.output, "parsed_output.txt")
    write_scenario_output(
        parsed_output_file, time_scale, scenarios, 
        args.mutation_rate, args.genome_length, args.generation_time
    )


if __name__ == "__main__":
    main()
