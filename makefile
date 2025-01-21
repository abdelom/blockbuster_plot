# Répertoires pour les sources, les objets et les exécutables
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# Noms des exécutables
TARGET1 = $(BIN_DIR)/blockbuster_grid_main
TARGET2 = $(BIN_DIR)/blockbuster_main

# Fichiers source (dans le répertoire src)
SRCS1 = $(SRC_DIR)/blockbuster_grid_main.c $(SRC_DIR)/blockbuster_grid.c
SRCS2 = $(SRC_DIR)/blockbuster_main.c $(SRC_DIR)/blockbuster.c $(SRC_DIR)/blockbuster_grid.c

# Fichiers objets générés (stockés dans le répertoire obj)
OBJS1 = $(OBJ_DIR)/blockbuster_grid_main.o $(OBJ_DIR)/blockbuster_grid.o
OBJS2 = $(OBJ_DIR)/blockbuster_main.o $(OBJ_DIR)/blockbuster.o $(OBJ_DIR)/blockbuster_grid.o

# Compilateur
CC = gcc

# Options de compilation
CFLAGS = -g -fopenmp

# Bibliothèques à lier
LIBS = -lm -lmpfr -llapacke -lopenblas

# Règle par défaut (compilation des deux exécutables)
all: $(TARGET1) $(TARGET2)

# Règle pour créer le premier exécutable (blockbuster_grid_main)
$(TARGET1): $(OBJS1)
	@mkdir -p $(BIN_DIR)  # Crée le répertoire bin si nécessaire
	$(CC) $(CFLAGS) -o $(TARGET1) $(OBJS1) $(LIBS)

# Règle pour créer le second exécutable (blockbuster_main)
$(TARGET2): $(OBJS2)
	@mkdir -p $(BIN_DIR)  # Crée le répertoire bin si nécessaire
	$(CC) $(CFLAGS) -o $(TARGET2) $(OBJS2) $(LIBS)

# Règle pour compiler les fichiers .c en fichiers .o (pour blockbuster_grid_main)
$(OBJ_DIR)/blockbuster_grid_main.o: $(SRC_DIR)/blockbuster_grid_main.c $(SRC_DIR)/blockbuster_grid.h
	@mkdir -p $(OBJ_DIR)  # Crée le répertoire obj si nécessaire
	$(CC) $(CFLAGS) -c $< -o $@

# Règle pour compiler blockbuster_main.c en objet
$(OBJ_DIR)/blockbuster_main.o: $(SRC_DIR)/blockbuster_main.c $(SRC_DIR)/blockbuster.h $(SRC_DIR)/blockbuster_grid.h
	@mkdir -p $(OBJ_DIR)  # Crée le répertoire obj si nécessaire
	$(CC) $(CFLAGS) -c $< -o $@

# Règle pour compiler blockbuster_grid.c en objet
$(OBJ_DIR)/blockbuster_grid.o: $(SRC_DIR)/blockbuster_grid.c $(SRC_DIR)/blockbuster_grid.h
	@mkdir -p $(OBJ_DIR)  # Crée le répertoire obj si nécessaire
	$(CC) $(CFLAGS) -c $< -o $@

# Règle pour compiler blockbuster.c en objet
$(OBJ_DIR)/blockbuster.o: $(SRC_DIR)/blockbuster.c $(SRC_DIR)/blockbuster.h
	@mkdir -p $(OBJ_DIR)  # Crée le répertoire obj si nécessaire
	$(CC) $(CFLAGS) -c $< -o $@

# Nettoyage des fichiers objets et des exécutables
clean:
	rm -f $(OBJS1) $(OBJS2) $(TARGET1) $(TARGET2)

# Forcer la recompilation
.PHONY: clean all
