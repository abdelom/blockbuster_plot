# Répertoires pour les sources, les objets et les exécutables
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# Noms des exécutables
TARGET1 = $(BIN_DIR)/blockbuster_grid_main
TARGET2 = $(BIN_DIR)/blockbuster_main
TARGET3 = $(BIN_DIR)/blockbuster_simulator

# Fichiers source
SRCS = $(wildcard $(SRC_DIR)/*.c)

# Conversion automatique en objets (src/foo.c → obj/foo.o)
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

# Dépendances spécifiques : chaque exécutable n’utilise pas tous les objets
OBJS1 = $(OBJ_DIR)/blockbuster_grid_main.o $(OBJ_DIR)/blockbuster_grid.o $(OBJ_DIR)/linear_regression.o $(OBJ_DIR)/linear_regression_f.o
OBJS2 = $(OBJ_DIR)/blockbuster_main.o $(OBJ_DIR)/blockbuster.o $(OBJ_DIR)/blockbuster_grid.o $(OBJ_DIR)/sfs.o $(OBJ_DIR)/linear_regression.o $(OBJ_DIR)/linear_regression_f.o
OBJS3 = $(OBJ_DIR)/blockbuster_simulator.o $(OBJ_DIR)/blockbuster_grid.o $(OBJ_DIR)/linear_regression.o $(OBJ_DIR)/linear_regression_f.o

# Compilateur et options
CC = gcc
CFLAGS = -g -fopenmp
LIBS = -lm -llapacke -lopenblas

# Règle par défaut
all: $(TARGET1) $(TARGET2) $(TARGET3)

# Règles de liens
$(TARGET1): $(OBJS1)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

$(TARGET2): $(OBJS2)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

$(TARGET3): $(OBJS3)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

# Règle générique de compilation (src/xxx.c → obj/xxx.o)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Nettoyage
clean:
	rm -f $(OBJ_DIR)/*.o $(TARGET1) $(TARGET2) $(TARGET3)

.PHONY: clean all
