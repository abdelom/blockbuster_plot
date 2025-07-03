import matplotlib.pyplot as plt

# Liste des valeurs fournies
values = [
    2461209, 1202462, 683586, 426382, 278390, 189684, 132721, 99468, 74577, 60957,
    50781, 43654, 40298, 35491, 33441, 32357, 30743, 31962, 33505, 36444,
    37754, 45709, 51765, 64169, 116148
]

# Multiplier chaque valeur par son index correspondant (en commençant par 1)
indexed_values = [value * (index + 1) for index, value in enumerate(values)]

# Tracer les résultats de la multiplication
plt.figure(figsize=(12, 6))
plt.plot(indexed_values, marker='o', linestyle='-')
plt.title('Graphique des valeurs multipliées par leur index')
plt.xlabel('Index')
plt.ylabel('Valeur multipliée par l\'index')
plt.grid(True)
plt.show()