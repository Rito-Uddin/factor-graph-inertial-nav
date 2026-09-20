# Fondements Théoriques : Variétés Riemanniennes, Groupes de Lie et Graphes de Facteurs

## 1. Introduction et Vue d'Ensemble
L'estimation d'état en robotique mobile (drones 6-DoF, navigation GNSS-denied) requiert la modélisation de grandeurs géométriques non linéaires, en particulier l'orientation 3D. 

L'approche classique repose sur le filtrage récursif de Kalman dans des espaces vectoriels $\mathbb{R}^n$. Cependant, l'utilisation de représentations minimales (angles d'Euler) introduit des singularités cinématiques (*gimbal lock*), tandis que les représentations redondantes (quaternions sur $S^3$) introduisent des contraintes de double recouvrement ($q \equiv -q$). 

La théorie de Lie permet de contourner ces limitations en traitant l'espace des poses comme une variété différentielle riemannienne, dotée d'un calcul différentiel rigoureux sans singularités.

---

## 2. Motivation : Pourquoi la Théorie de Lie en Robotique ?

### Fondement de la théorie
Comme le soulignent Joan Solà et al. dans *A micro Lie theory for state estimation in robotics* :
> *« Relying on the Lie theory, we are able to construct a rigorous calculus corpus to handle uncertainties, derivatives, and integrals with precision and ease. [...] Amazingly, the group $\mathcal{G}$ is almost completely determined by $\mathfrak{g}$ and its Lie bracket. Thus for many purposes one can replace $\mathcal{G}$ with $\mathfrak{g}$. Since $\mathcal{G}$ is a complicated nonlinear object and $\mathfrak{g}$ is just a vector space, it is usually vastly simpler to work with $\mathfrak{g}$. »*

En robotique, il est d'usage pragmatique d'exploiter directement l'espace vectoriel tangent $\mathbb{R}^n$ isomorphe à l'algèbre $\mathfrak{g}$ (via les opérateurs $[\cdot]_\times$ et $(\cdot)^\vee$), ce qui permet de définir des dérivées, covariances et incertitudes sans manipuler lourdement les crochets de Lie continus.

### Propriété d'homogénéité des groupes de Lie
Une variété différentiable (*smooth manifold*) est un espace topologique qui ressemble localement à un espace euclidien. Sur un groupe de Lie $\mathcal{G}$, la structure est homogène : le voisinage de n'importe quel élément ressemble au voisinage de l'élément neutre (l'identité $\mathcal{E}$). 

Par conséquent, tous les espaces tangents en tout point sont isomorphes à l'espace tangent à l'identité, qui définit l'algèbre de Lie $\mathfrak{g}$. Le repère d'origine est identifié à l'identité du groupe, et tout point de la variété représente un repère local particulier.

### Comparaison des paradigmes d'estimation

| Méthode | Gestion du passé | Linéarisation | Comportement blackout |
| :--- | :--- | :--- | :--- |
| **EKF standard** | Oubli récursif | Précoce (figée à $t-1$) | Divergence, saut discontinu |
| **UKF** | Oubli récursif | Points sigma dans $\mathbb{R}^n$ | Rupture géométrique sur $SO(3)$ |
| **Graphe (iSAM2)** | Lissage de la trajectoire | Rétroactive fluide | Re-convergence globale |

### Groupes de Lie matriciels usuels en robotique

| Variété / Groupe | Élément $X$ | Espace tangent $\mathbb{R}^n$ | Dimension | Description |
| :--- | :---: | :---: | :---: | :--- |
| $SO(3)$ | $R \in \mathbb{R}^{3 \times 3}$ | $\omega \in \mathbb{R}^3$ | 3 | Rotations 3D |
| $SE(3)$ | $T = \begin{bmatrix} R & p \\ 0 & 1 \end{bmatrix}$ | $\xi = \begin{bmatrix} \rho \\ \theta \end{bmatrix} \in \mathbb{R}^6$ | 6 | Poses rigides 3D (rotation + translation) |

---

## 3. Définitions et Notations sur $SO(3)$

### Le Groupe Spécial Orthogonal
Le groupe $SO(3)$ modélise l'attitude d'un corps rigide sans aucune singularité :
$$SO(3) = \left\{ R \in \mathbb{R}^{3 \times 3} \;\middle\vert{}\; R^T R = I_3, \; \det(R) = +1 \right\}$$

$SO(3)$ est un groupe pour la multiplication matricielle mais n'est **pas** un espace vectoriel ($R_1 + R_2 \notin SO(3)$). L'inverse d'une rotation est son adjointe : $R^{-1} = R^T$.

### Algèbre associée $\mathfrak{so}(3)$ et Opérateur Hat
L'algèbre de Lie $\mathfrak{so}(3)$ est l'espace vectoriel tangent à l'identité, constitué des matrices réelles antisymétriques d'ordre 3 :
$$\mathfrak{so}(3) = \left\{ \Omega \in \mathbb{R}^{3 \times 3} \;\middle\vert{}\; \Omega^T = -\Omega \right\}$$

L'opérateur *hat* $[\cdot]_\times : \mathbb{R}^3 \to \mathfrak{so}(3)$ associe un vecteur de vitesse angulaire à sa matrice antisymétrique :
$$\omega = \begin{bmatrix} \omega_1 \\ \omega_2 \\ \omega_3 \end{bmatrix} \implies [\omega]_\times = \begin{bmatrix} 0 & -\omega_3 & \omega_2 \\ \omega_3 & 0 & -\omega_1 \\ -\omega_2 & \omega_1 & 0 \end{bmatrix}$$

L'opération réciproque est l'opérateur *vee* noté $(\cdot)^\vee : \mathfrak{so}(3) \to \mathbb{R}^3$ tel que $([\omega]_\times)^\vee = \omega$.

### Cartes Exponentielle et Logarithmique
Le passage rigoureux entre l'algèbre $\mathfrak{so}(3)$ (locale, linéaire) et le groupe $SO(3)$ (global, non linéaire) s'effectue via :

* **Carte exponentielle ($\exp$)** : Donnée analytiquement par la formule de Rodrigues pour $\theta = \Vert{}\phi\Vert{}_2$ et l'axe unitaire $u = \phi / \theta$ :
  $$\exp([\phi]_\times) = I_3 + \frac{\sin \theta}{\theta}[\phi]_\times + \frac{1 - \cos \theta}{\theta^2}[\phi]_\times^2$$

* **Carte logarithmique ($\log$)** : Permet d'extraire le vecteur axe-angle depuis une matrice $R \in SO(3)$ :
  $$\theta = \arccos\left(\frac{\mathrm{tr}(R) - 1}{2}\right), \quad [\phi]_\times = \frac{\theta}{2\sin \theta}(R - R^T)$$

### Opérateurs de perturbation locale $\boxplus$ et $\boxminus$
Dans l'optimiseur, une correction locale $\delta\phi \in \mathbb{R}^3$ est injectée sur une rotation estimée $R$ via la rétraction locale à droite (repère corps/local) :
$$R \boxplus \delta\phi = R \cdot \exp([\delta\phi]_\times)$$

L'écart tangent (résidu géodésique) entre deux rotations s'exprime par :
$$R_1 \boxminus R_2 = \log(R_2^T R_1)^\vee \in \mathbb{R}^3$$

---

## 4. Développements et Préintégration IMU

### Cinématique IMU en temps continu
Les mesures fournies par les accéléromètres $\tilde{a}(t)$ et gyromètres $\tilde{\omega}(t)$ s'écrivent :
$$\tilde{\omega}(t) = \omega(t) + b_g(t) + \eta_g(t)$$
$$\tilde{a}(t) = R^T(t) \left( a_{\text{world}}(t) - g \right) + b_a(t) + \eta_a(t)$$

où $\eta_a, \eta_g$ sont des bruits blancs gaussiens et $b_a, b_g$ suivent des marches aléatoires dérivantes (*Brownian motion*).

### Principe de la préintégration (Forster et al.)
Pour éviter de réintégrer l'ensemble des mesures inertielles à haute cadence lors des itérations du solveur d'optimisation, les incréments de mouvement sont formulés dans le repère local de l'instant $i$ :
$$\Delta R_{ij} = \prod_{k=i}^{j-1} \exp\left( [\tilde{\omega}_k - \bar{b}_g^i]_\times \Delta t \right)$$
$$\Delta v_{ij} = \sum_{k=i}^{j-1} \Delta R_{ik} (\tilde{a}_k - \bar{b}_a^i) \Delta t$$
$$\Delta p_{ij} = \sum_{k=i}^{j-1} \left( \Delta v_{ik}\Delta t + \frac{1}{2}\Delta R_{ik}(\tilde{a}_k - \bar{b}_a^i)\Delta t^2 \right)$$

### Correction linéaire par Jacobiennes de biais
Lors de la mise à jour des biais $\delta b = b - \bar{b}$ par l'optimiseur, une expansion de Taylor au premier ordre actualise les deltas sans réintégration numérique :
$$\Delta R_{ij}(b_g) \approx \Delta R_{ij}(\bar{b}_g) \cdot \exp\left( \left[ \frac{\partial \Delta R_{ij}}{\partial b_g} \delta b_g \right]_\times \right)$$

---

## 5. Notes d'Implémentation sous GTSAM
* L'attitude est représentée par la classe `gtsam.Rot3`, avec `Rot3.Expmap` et `Rot3.Logmap`.
* Les états successifs sont indexés par les symboles shorthand : `X(k)` pour la pose $SE(3)$, `V(k)` pour la vitesse $\mathbb{R}^3$, et `B(k)` pour les biais inertiels.
* La préintégration IMU est confiée au composant `gtsam.PreintegratedCombinedMeasurements`.
* L'optimisation incrémentale temps-réel est résolue par `gtsam.ISAM2` sur la structure d'arbre de Bayes (*Bayes Tree*).

---

## 6. Références
1. J. Solà, J. Deray, and D. Atchuthan, *« A micro Lie theory for state estimation in robotics »*, arXiv:1812.01537, 2018.
2. C. Forster, L. Carlone, F. Dellaert, and D. Scaramuzza, *« IMU Preintegration on Manifold for High-Frame-Rate GNSS-denied Robotics »*, IEEE Transactions on Robotics, vol. 33, no. 1, pp. 249–265, 2017.
3. F. Dellaert and M. Kaess, *« Factor Graphs for Robot Perception »*, Foundations and Trends in Robotics, vol. 6, no. 1-2, pp. 1–139, 2017.