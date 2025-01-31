

#!/bin/bash

# Liste des commandes à exécuter
commands=(
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Athene_cunicularia.fs -p ../inf5/Athene_cunicularia/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 1088737051 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Corvus_brachyrhynchos.fs -p ../inf5/Corvus_brachyrhynchos/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 378603033 -g 3 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Corvus_corone_corone.fs -p ../inf5/Corvus_corone_corone/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 472795516 -g 2 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Corvus_dauuricus.fs -p ../inf5/Corvus_dauuricus/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 976409105 -g 2 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Corvus_frugilegus.fs -p ../inf5/Corvus_frugilegus/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 893602510 -g 2 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Coturnix_japonica.fs -p ../inf5/Coturnix_japonica/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 483720035 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Egretta_garzetta.fs -p ../inf5/Egretta_garzetta/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 491257580 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Eopsaltria_australis.fs -p ../inf5/Eopsaltria_australis/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 648764071 -g 6.3 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Falco_peregrinus.fs -p ../inf5/Falco_peregrinus/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 571549750 -g 2 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_albicollis.fs -p ../inf5/Ficedula_albicollis/ -o 0 -c 8 -l 1e-4 -b 6 -m 4.60E-09 -L 1004991485 -g 1 "
)


# Utilisation de parallel pour exécuter les commandes en parallèle, avec une limite de 8 jobs simultanés
parallel -j 8 ::: "${commands[@]}"



python ../rapports.py --sfs_file ../SFS/SFS_Athene_cunicularia.fs --aic_plot ../inf5/Athene_cunicularia/aic_plot.png --likelihood_plot ../inf5/Athene_cunicularia/log_likelihood_plot.png --output ../rapports/Athene_cunicularia/

python ../rapports.py --sfs_file ../SFS/SFS_Corvus_brachyrhynchos.fs --aic_plot ../inf5/Corvus_brachyrhynchos/aic_plot.png --likelihood_plot ../inf5/Corvus_brachyrhynchos/log_likelihood_plot.png --output ../rapports/Corvus_brachyrhynchos/

python ../rapports.py --sfs_file ../SFS/SFS_Corvus_corone_corone.fs --aic_plot ../inf5/Corvus_corone_corone/aic_plot.png --likelihood_plot ../inf5/Corvus_corone_corone/log_likelihood_plot.png --output ../rapports/Corvus_corone_corone/

python ../rapports.py --sfs_file ../SFS/SFS_Corvus_dauuricus.fs --aic_plot ../inf5/Corvus_dauuricus/aic_plot.png --likelihood_plot ../inf5/Corvus_dauuricus/log_likelihood_plot.png --output ../rapports/Corvus_dauuricus/

python ../rapports.py --sfs_file ../SFS/SFS_Corvus_frugilegus.fs --aic_plot ../inf5/Corvus_frugilegus/aic_plot.png --likelihood_plot ../inf5/Corvus_frugilegus/log_likelihood_plot.png --output ../rapports/Corvus_frugilegus/

python ../rapports.py --sfs_file ../SFS/SFS_Coturnix_japonica.fs --aic_plot ../inf5/Coturnix_japonica/aic_plot.png --likelihood_plot ../inf5/Coturnix_japonica/log_likelihood_plot.png --output ../rapports/Coturnix_japonica/

python ../rapports.py --sfs_file ../SFS/SFS_Eopsaltria_australis.fs --aic_plot ../inf5/Eopsaltria_australis/aic_plot.png --likelihood_plot ../inf5/Eopsaltria_australis/log_likelihood_plot.png --output ../rapports/Eopsaltria_australis/

python ../rapports.py --sfs_file ../SFS/SFS_Falco_peregrinus.fs --aic_plot ../inf5/Falco_peregrinus/aic_plot.png --likelihood_plot ../inf5/Falco_peregrinus/log_likelihood_plot.png --output ../rapports/Falco_peregrinus/

python ../rapports.py --sfs_file ../SFS/SFS_Ficedula_albicollis.fs --aic_plot ../inf5/Ficedula_albicollis/aic_plot.png --likelihood_plot ../inf5/Ficedula_albicollis/log_likelihood_plot.png --output ../rapports/Ficedula_albicollis/


# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Coturnix_japonica.fs -p ../inft/Coturnix_japonica/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 483720035 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Egretta_garzetta.fs -p ../inft/Egretta_garzetta/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 491257580 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Eopsaltria_australis.fs -p ../inft/Eopsaltria_australis/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 648764071 -g 6.3 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Falco_peregrinus.fs -p ../inft/Falco_peregrinus/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 571549750 -g 2 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_albicollis.fs -p ../inft/Ficedula_albicollis/ -o 0 -c 8 -l 1e-5 -b 6 -m 4.60E-09 -L 1004991485 -g 1  "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_hypoleuca.fs -p ../inft/Ficedula_hypoleuca/ -o 0 -c 8 -l 1e-5 -b 6 -m 5.00E-09 -L 666193910 -g 1.5 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_semitorquata.fs -p ../inft/Ficedula_semitorquata/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 624431834 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_speculigera.fs -p ../inft/Ficedula_speculigera/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 558313357 -g 1.5 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Parus_major.fs -p ../inft/Parus_major/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 777237111 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Parus_major_eur.fs -p ../inft/Parus_major_eur/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 777237111 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Passer_domesticus.fs -p ../inft/Passer_domesticus/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 591106486 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Passer_hispaniolensis.fs -p ../inft/Passer_hispaniolensis/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 587601377 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Passer_italiae.fs -p ../inft/Passer_italiae/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 674510400 -g 1"
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Passer_montanus.fs -p ../inft/Passer_montanus/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 506071734 -g 1"
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Phylloscopus_collybita_abietinus.fs -p ../inft/Phylloscopus_collybita_abietinus/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 894720137 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Phylloscopus_collybita_tristis.fs -p ../inft/Phylloscopus_collybita_tristis/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 944155683 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Phylloscopus_trochilus_acredula.fs -p ../inft/Phylloscopus_trochilus_acredula/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 546825861 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Poephila_acuticauda.fs -p ../inft/Poephila_acuticauda/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 555992342 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Sylvia_atricapilla.fs -p ../inft/Sylvia_atricapilla/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 661731274 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Sylvia_borin.fs -p ../inft/Sylvia_borin/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 487726329 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Taeniopygia_guttata.fs -p ../inft/Taeniopygia_guttata/ -o 0 -c 8 -l 1e-4 -b 6 -m 6.63E-09 -L 625640347 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Hirundo_rustica_mnhn_contemp.fs -p ../inft/Hirundo_rustica_mnhn_contemp/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 988409861 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Probosciger_aterrimus.fs -p ../inft/Probosciger_aterrimus/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 313079153 -g 3 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Muscicapa_striata_striata.fs -p ../inft/Muscicapa_striata_striata/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 466388329 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Picus_viridis.fs -p ../inft/Picus_viridis/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 712348313 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Nipponia_nippon.fs -p ../inft/Nipponia_nippon/ -o 0 -c 8 -l 1e-4 -b 6 -m 5.00E-09 -L 328377865 -g 4 "