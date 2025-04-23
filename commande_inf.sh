

#!/bin/bash

# Liste des commandes à exécuter
commands=(
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Athene_cunicularia.fs -p ../infbird/Athene_cunicularia/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 1088737051 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Corvus_brachyrhynchos.fs -p ../infbird/Corvus_brachyrhynchos/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 378603033 -g 3 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Corvus_corone_corone.fs -p ../infbird/Corvus_corone_corone/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 472795516 -g 2 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Corvus_dauuricus.fs -p ../infbird/Corvus_dauuricus/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 976409105 -g 2 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Corvus_frugilegus.fs -p ../infbird/Corvus_frugilegus/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 893602510 -g 2 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Coturnix_japonica.fs -p ../infbird/Coturnix_japonica/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 483720035 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Egretta_garzetta.fs -p ../infbird/Egretta_garzetta/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 491257580 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Eopsaltria_australis.fs -p ../infbird/Eopsaltria_australis/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 648764071 -g 6.3 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Falco_peregrinus.fs -p ../infbird/Falco_peregrinus/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 571549750 -g 2 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_albicollis.fs -p ../infbird/Ficedula_albicollis/ -o 0 -c 6 -n 40 -b 3     -m 4.60E-09 -L 1004991485 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_hypoleuca.fs -p ../infbird/Ficedula_hypoleuca/ -o 0 -c 6 -n 40 -b 3   -l 1e-5  -m 5.00E-09 -L 666193910 -g 1.5 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_semitorquata.fs -p ../infbird/Ficedula_semitorquata/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 624431834 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_speculigera.fs -p ../infbird/Ficedula_speculigera/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 558313357 -g 1.5 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Parus_major.fs -p ../infbird/Parus_major/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 777237111 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Parus_major_eur.fs -p ../infbird/Parus_major_eur/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 777237111 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Passer_domesticus.fs -p ../infbird/Passer_domesticus/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 591106486 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Passer_hispaniolensis.fs -p ../infbird/Passer_hispaniolensis/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 587601377 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Passer_italiae.fs -p ../infbird/Passer_italiae/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 674510400 -g 1"
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Passer_montanus.fs -p ../infbird/Passer_montanus/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 506071734 -g 1"
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Phylloscopus_collybita_abietinus.fs -p ../infbird/Phylloscopus_collybita_abietinus/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 894720137 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Phylloscopus_collybita_tristis.fs -p ../infbird/Phylloscopus_collybita_tristis/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 944155683 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Phylloscopus_trochilus_acredula.fs -p ../infbird/Phylloscopus_trochilus_acredula/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 546825861 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Poephila_acuticauda.fs -p ../infbird/Poephila_acuticauda/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 555992342 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Sylvia_atricapilla.fs -p ../infbird/Sylvia_atricapilla/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 661731274 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Sylvia_borin.fs -p ../infbird/Sylvia_borin/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 487726329 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Taeniopygia_guttata.fs -p ../infbird/Taeniopygia_guttata/ -o 0 -c 6 -n 40 -b 3     -m 6.63E-09 -L 625640347 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Hirundo_rustica_mnhn_contemp.fs -p ../infbird/Hirundo_rustica_mnhn_contemp/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 988409861 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Probosciger_aterrimus.fs -p ../infbird/Probosciger_aterrimus/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 313079153 -g 3 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Muscicapa_striata_striata.fs -p ../infbird/Muscicapa_striata_striata/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 466388329 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Picus_viridis.fs -p ../infbird/Picus_viridis/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 712348313 -g 1 "
"bash blockbuster_plot.sh --sfs ../SFS/SFS_Nipponia_nippon.fs -p ../infbird/Nipponia_nippon/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 328377865 -g 4 "
)


# Utilisation de parallel pour exécuter les commandes en parallèle, avec une limite de 8 jobs simultanés
parallel -j 14 ::: "${commands[@]}"



python ../rapports/rapports.py --sfs_file ../SFS/SFS_Athene_cunicularia.fs --aic_plot ../inf5/Athene_cunicularia/aic_plot.png --likelihood_plot ../inf5/Athene_cunicularia/log_likelihood_plot.png --output ../rapports/Athene_cunicularia/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Corvus_brachyrhynchos.fs --aic_plot ../inf5/Corvus_brachyrhynchos/aic_plot.png --likelihood_plot ../inf5/Corvus_brachyrhynchos/log_likelihood_plot.png --output ../rapports/Corvus_brachyrhynchos/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Corvus_corone_corone.fs --aic_plot ../inf5/Corvus_corone_corone/aic_plot.png --likelihood_plot ../inf5/Corvus_corone_corone/log_likelihood_plot.png --output ../rapports/Corvus_corone_corone/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Corvus_dauuricus.fs --aic_plot ../inf5/Corvus_dauuricus/aic_plot.png --likelihood_plot ../inf5/Corvus_dauuricus/log_likelihood_plot.png --output ../rapports/Corvus_dauuricus/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Corvus_frugilegus.fs --aic_plot ../inf5/Corvus_frugilegus/aic_plot.png --likelihood_plot ../inf5/Corvus_frugilegus/log_likelihood_plot.png --output ../rapports/Corvus_frugilegus/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Coturnix_japonica.fs --aic_plot ../inf5/Coturnix_japonica/aic_plot.png --likelihood_plot ../inf5/Coturnix_japonica/log_likelihood_plot.png --output ../rapports/Coturnix_japonica/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Egretta_garzetta.fs --aic_plot ../inf5/Egretta_garzetta/aic_plot.png --likelihood_plot ../inf5/Egretta_garzetta/log_likelihood_plot.png --output ../rapports/Egretta_garzetta/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Eopsaltria_australis.fs --aic_plot ../inf5/Eopsaltria_australis/aic_plot.png --likelihood_plot ../inf5/Eopsaltria_australis/log_likelihood_plot.png --output ../rapports/Eopsaltria_australis/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Falco_peregrinus.fs --aic_plot ../inf5/Falco_peregrinus/aic_plot.png --likelihood_plot ../inf5/Falco_peregrinus/log_likelihood_plot.png --output ../rapports/Falco_peregrinus/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Ficedula_albicollis.fs --aic_plot ../inf5/Ficedula_albicollis/aic_plot.png --likelihood_plot ../inf5/Ficedula_albicollis/log_likelihood_plot.png --output ../rapports/Ficedula_albicollis/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Ficedula_hypoleuca.fs --aic_plot ../inf5/Ficedula_hypoleuca/aic_plot.png --likelihood_plot ../inf5/Ficedula_hypoleuca/log_likelihood_plot.png --output ../rapports/Ficedula_hypoleuca/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Ficedula_semitorquata.fs --aic_plot ../inf5/Ficedula_semitorquata/aic_plot.png --likelihood_plot ../inf5/Ficedula_semitorquata/log_likelihood_plot.png --output ../rapports/Ficedula_semitorquata/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Ficedula_speculigera.fs --aic_plot ../inf5/Ficedula_speculigera/aic_plot.png --likelihood_plot ../inf5/Ficedula_speculigera/log_likelihood_plot.png --output ../rapports/Ficedula_speculigera/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Parus_major.fs --aic_plot ../inf5/Parus_major/aic_plot.png --likelihood_plot ../inf5/Parus_major/log_likelihood_plot.png --output ../rapports/Parus_major/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Parus_major_eur.fs --aic_plot ../inf5/Parus_major_eur/aic_plot.png --likelihood_plot ../inf5/Parus_major_eur/log_likelihood_plot.png --output ../rapports/Parus_major_eur/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Passer_domesticus.fs --aic_plot ../inf5/Passer_domesticus/aic_plot.png --likelihood_plot ../inf5/Passer_domesticus/log_likelihood_plot.png --output ../rapports/Passer_domesticus/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Passer_hispaniolensis.fs --aic_plot ../inf5/Passer_hispaniolensis/aic_plot.png --likelihood_plot ../inf5/Passer_hispaniolensis/log_likelihood_plot.png --output ../rapports/Passer_hispaniolensis/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Passer_italiae.fs --aic_plot ../inf5/Passer_italiae/aic_plot.png --likelihood_plot ../inf5/Passer_italiae/log_likelihood_plot.png --output ../rapports/Passer_italiae/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Passer_montanus.fs --aic_plot ../inf5/Passer_montanus/aic_plot.png --likelihood_plot ../inf5/Passer_montanus/log_likelihood_plot.png --output ../rapports/Passer_montanus/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Phylloscopus_collybita_abietinus.fs --aic_plot ../inf5/Phylloscopus_collybita_abietinus/aic_plot.png --likelihood_plot ../inf5/Phylloscopus_collybita_abietinus/log_likelihood_plot.png --output ../rapports/Phylloscopus_collybita_abietinus/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Phylloscopus_collybita_tristis.fs --aic_plot ../inf5/Phylloscopus_collybita_tristis/aic_plot.png --likelihood_plot ../inf5/Phylloscopus_collybita_tristis/log_likelihood_plot.png --output ../rapports/Phylloscopus_collybita_tristis/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Phylloscopus_trochilus_acredula.fs --aic_plot ../inf5/Phylloscopus_trochilus_acredula/aic_plot.png --likelihood_plot ../inf5/Phylloscopus_trochilus_acredula/log_likelihood_plot.png --output ../rapports/Phylloscopus_trochilus_acredula/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Poephila_acuticauda.fs --aic_plot ../inf5/Poephila_acuticauda/aic_plot.png --likelihood_plot ../inf5/Poephila_acuticauda/log_likelihood_plot.png --output ../rapports/Poephila_acuticauda/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Sylvia_atricapilla.fs --aic_plot ../inf5/Sylvia_atricapilla/aic_plot.png --likelihood_plot ../inf5/Sylvia_atricapilla/log_likelihood_plot.png --output ../rapports/Sylvia_atricapilla/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Sylvia_borin.fs --aic_plot ../inf5/Sylvia_borin/aic_plot.png --likelihood_plot ../inf5/Sylvia_borin/log_likelihood_plot.png --output ../rapports/Sylvia_borin/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Taeniopygia_guttata.fs --aic_plot ../inf5/Taeniopygia_guttata/aic_plot.png --likelihood_plot ../inf5/Taeniopygia_guttata/log_likelihood_plot.png --output ../rapports/Taeniopygia_guttata/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Hirundo_rustica_mnhn_contemp.fs --aic_plot ../inf5/Hirundo_rustica_mnhn_contemp/aic_plot.png --likelihood_plot ../inf5/Hirundo_rustica_mnhn_contemp/log_likelihood_plot.png --output ../rapports/Hirundo_rustica_mnhn_contemp/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Probosciger_aterrimus.fs --aic_plot ../inf5/Probosciger_aterrimus/aic_plot.png --likelihood_plot ../inf5/Probosciger_aterrimus/log_likelihood_plot.png --output ../rapports/Probosciger_aterrimus/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Muscicapa_striata_striata.fs --aic_plot ../inf5/Muscicapa_striata_striata/aic_plot.png --likelihood_plot ../inf5/Muscicapa_striata_striata/log_likelihood_plot.png --output ../rapports/Muscicapa_striata_striata/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Picus_viridis.fs --aic_plot ../inf5/Picus_viridis/aic_plot.png --likelihood_plot ../inf5/Picus_viridis/log_likelihood_plot.png --output ../rapports/Picus_viridis/

python ../rapports/rapports.py --sfs_file ../SFS/SFS_Nipponia_nippon.fs --aic_plot ../inf5/Nipponia_nippon/aic_plot.png --likelihood_plot ../inf5/Nipponia_nippon/log_likelihood_plot.png --output ../rapports/Nipponia_nippon/

# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Coturnix_japonica.fs -p ../infbird/Coturnix_japonica/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 483720035 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Egretta_garzetta.fs -p ../infbird/Egretta_garzetta/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 491257580 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Eopsaltria_australis.fs -p ../infbird/Eopsaltria_australis/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 648764071 -g 6.3 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Falco_peregrinus.fs -p ../infbird/Falco_peregrinus/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 571549750 -g 2 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_albicollis.fs -p ../infbird/Ficedula_albicollis/ -o 0 -c 6 -n 40 -b 3   -l 1e-5  -m 4.60E-09 -L 1004991485 -g 1  "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_hypoleuca.fs -p ../infbird/Ficedula_hypoleuca/ -o 0 -c 6 -n 40 -b 3   -l 1e-5  -m 5.00E-09 -L 666193910 -g 1.5 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_semitorquata.fs -p ../infbird/Ficedula_semitorquata/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 624431834 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Ficedula_speculigera.fs -p ../infbird/Ficedula_speculigera/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 558313357 -g 1.5 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Parus_major.fs -p ../infbird/Parus_major/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 777237111 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Parus_major_eur.fs -p ../infbird/Parus_major_eur/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 777237111 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Passer_domesticus.fs -p ../infbird/Passer_domesticus/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 591106486 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Passer_hispaniolensis.fs -p ../infbird/Passer_hispaniolensis/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 587601377 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Passer_italiae.fs -p ../infbird/Passer_italiae/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 674510400 -g 1"
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Passer_montanus.fs -p ../infbird/Passer_montanus/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 506071734 -g 1"
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Phylloscopus_collybita_abietinus.fs -p ../infbird/Phylloscopus_collybita_abietinus/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 894720137 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Phylloscopus_collybita_tristis.fs -p ../infbird/Phylloscopus_collybita_tristis/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 944155683 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Phylloscopus_trochilus_acredula.fs -p ../infbird/Phylloscopus_trochilus_acredula/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 546825861 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Poephila_acuticauda.fs -p ../infbird/Poephila_acuticauda/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 555992342 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Sylvia_atricapilla.fs -p ../infbird/Sylvia_atricapilla/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 661731274 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Sylvia_borin.fs -p ../infbird/Sylvia_borin/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 487726329 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Taeniopygia_guttata.fs -p ../infbird/Taeniopygia_guttata/ -o 0 -c 6 -n 40 -b 3     -m 6.63E-09 -L 625640347 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Hirundo_rustica_mnhn_contemp.fs -p ../infbird/Hirundo_rustica_mnhn_contemp/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 988409861 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Probosciger_aterrimus.fs -p ../infbird/Probosciger_aterrimus/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 313079153 -g 3 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Muscicapa_striata_striata.fs -p ../infbird/Muscicapa_striata_striata/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 466388329 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Picus_viridis.fs -p ../infbird/Picus_viridis/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 712348313 -g 1 "
# "bash blockbuster_plot.sh --sfs ../SFS/SFS_Nipponia_nippon.fs -p ../infbird/Nipponia_nippon/ -o 0 -c 6 -n 40 -b 3     -m 5.00E-09 -L 328377865 -g 4 "