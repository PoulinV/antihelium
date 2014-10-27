cmake_minimum_required (VERSION 2.8)


#Declaration du projet

project(antiprotons)

#Declaration des fichiers sources

set(SOURCES_FILES nrutil.c GAUSSJ.c TRIDAG.c besselj0_next.c besselj1_next.c BESSEL_PRELIM.c CROSS_SECTIONS.c DIFFUSION_PROPAGATION.c SOLAR_MOD.c PROTON.c HELIUM.c PRIMARY_PBAR.c ANTI_PROTON.c MAIN.c)

#Recherche des librairies sur le systeme

find_library(M_LIB m)

#Declaration de l'executable

add_executable(executable ${SOURCES_FILES})

#Lien avec les libraries

target_link_libraries(executable ${M_LIB})