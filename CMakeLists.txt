cmake_minimum_required (VERSION 2.8)


#Declaration du projet

project(antiprotons)

#Declaration des fichiers sources

set(SOURCES_FILES sources/nrutil.c sources/GAUSSJ.c sources/TRIDAG.c sources/besselj0_next.c sources/besselj1_next.c sources/BESSEL_PRELIM.c sources/CROSS_SECTIONS.c sources/DIFFUSION_PROPAGATION.c sources/SOLAR_MOD.c sources/PROTON.c sources/HELIUM.c sources/PRIMARY_PBAR.c sources/ANTI_PROTON.c MAIN.c)

#Recherche des librairies sur le systeme

find_library(M_LIB m)

#Declaration de l'executable

add_executable(executable ${SOURCES_FILES})

#Lien avec les libraries

target_link_libraries(executable ${M_LIB})