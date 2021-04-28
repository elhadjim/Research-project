Pour les méthodes LBP, Eigenfaces, Fisherfaces et les réseaux de neurones d'Elhadji et AlexNet, il y a un script qui fait l'apprentissage et une fonction qui fait la reconnaissance/détection.
Pour générer les données d'apprentissage, il faut d'abord exécuter  les scripts suivants (dans le même dossier que celui où se trouve le fichier pir_gui.m) :
- LBP_EyeTrainingSetConstruction.m
- LBP_TrainingSetConstruction.m
- Eigenfaces_TrainingSetConstruction.m
- Fisherfaces_TrainingSetConstruction.m
- AlexNet_Training.m
- CNNElhadji_Training.m

Pour lancer l'interface graphique, il suffit ensuite d'exécuter le script suivant :
pir_gui.m

L'interface permet à la fois de prendre une image avec la caméra ou d'importer un fichier déjà existant sur l'ordinateur. On peut alors lancer la détection avec LBP ou le réseau de neurones.
Pour lancer la reconnaissance, il faut avoir lancé la détection au préalable si on a pris une photo avec la caméra. La reconnaissance s'effectuera alors sur les visages détectés lors de la détection.
Sinon, si on importe une image existante, on peut lancer la reconnaissance sans avoir lancé la détection et dans ce cas, la reconnaissance sera effectuée sur l'image entière (il faut donc que celle-ci soit bien cadrée, et si possible en format carré ou bien portrait 4:3).

Pour SIFT et SURF, les algorithmes nécessitent deux photos, il faut donc suivre les indications données par l'interface.

Pour SURF, la première photo (celle que l'on prend avec les boutons à gauche de l'interface) est celle qui doit contenir l'objet que l'on veut détecter dans une scène plus globale, tandis que la deuxième photo (celle que l'on prend avec les boutons qui apparaissent à droite de l'interface) est celle qui doit contenir seulement l'objet à détecter, en gros plan.