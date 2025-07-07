# Répertoires pour les sources, les objets et les exécutables
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# Noms des exécutables
TARGET1 = $(BIN_DIR)/blockbuster_grid_main
TARGET2 = $(BIN_DIR)/blockbuster_main
TARGET3 = $(BIN_DIR)/blockbuster_simulator
TARGET4 = $(BIN_DIR)/blockbuster_refine

# Fichiers source (dans le répertoire src)
SRCS1 = $(SRC_DIR)/blockbuster_grid_main.c $(SRC_DIR)/blockbuster_grid.c
SRCS2 = $(SRC_DIR)/blockbuster_main.c $(SRC_DIR)/blockbuster.c $(SRC_DIR)/blockbuster_grid.c
SRCS3 = $(SRC_DIR)/blockbuster_simulator.c $(SRC_DIR)/blockbuster_grid.c
SRCS4 = $(SRC_DIR)/blockbuster_refine.c $(SRC_DIR)/blockbuster.c $(SRC_DIR)/blockbuster_grid.c

# Fichiers objets générés (stockés dans le répertoire obj)
OBJS1 = $(OBJ_DIR)/blockbuster_grid_main.o $(OBJ_DIR)/blockbuster_grid.o
OBJS2 = $(OBJ_DIR)/blockbuster_main.o $(OBJ_DIR)/blockbuster.o $(OBJ_DIR)/blockbuster_grid.o
OBJS3 = $(OBJ_DIR)/blockbuster_simulator.o $(OBJ_DIR)/blockbuster_grid.o
OBJS4 = $(OBJ_DIR)/blockbuster_refine.o $(OBJ_DIR)/blockbuster.o $(OBJ_DIR)/blockbuster_grid.o

# Compilateur
CC = gcc

# Options de compilation
CFLAGS = -g -fopenmp

# Bibliothèques à lier
LIBS = -lm -llapacke -lopenblas

# Règle par défaut (compilation des quatre exécutables)
all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

# Règle pour créer le premier exécutable (blockbuster_grid_main)
$(TARGET1): $(OBJS1)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(TARGET1) $(OBJS1) $(LIBS)

# Règle pour créer le second exécutable (blockbuster_main)
$(TARGET2): $(OBJS2)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(TARGET2) $(OBJS2) $(LIBS)

# Règle pour créer le troisième exécutable (blockbuster_simulator)
$(TARGET3): $(OBJS3)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(TARGET3) $(OBJS3) $(LIBS)

# Règle pour créer le quatrième exécutable (blockbuster_refine)
$(TARGET4): $(OBJS4)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(TARGET4) $(OBJS4) $(LIBS)

# Règles de compilation des fichiers objets

$(OBJ_DIR)/blockbuster_grid_main.o: $(SRC_DIR)/blockbuster_grid_main.c $(SRC_DIR)/blockbuster_grid.h
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/blockbuster_main.o: $(SRC_DIR)/blockbuster_main.c $(SRC_DIR)/blockbuster.h $(SRC_DIR)/blockbuster_grid.h
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/blockbuster_grid.o: $(SRC_DIR)/blockbuster_grid.c $(SRC_DIR)/blockbuster_grid.h
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/blockbuster.o: $(SRC_DIR)/blockbuster.c $(SRC_DIR)/blockbuster.h
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/blockbuster_simulator.o: $(SRC_DIR)/blockbuster_simulator.c $(SRC_DIR)/blockbuster_grid.h
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/blockbuster_refine.o: $(SRC_DIR)/blockbuster_refine.c $(SRC_DIR)/blockbuster.h $(SRC_DIR)/blockbuster_grid.h
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Nettoyage des fichiers objets et des exécutables
clean:
	rm -f $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

# Forcer la recompilation
.PHONY: clean all
