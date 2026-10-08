"""Deterministic stem/loop diagrams of nested pairs, not molecular coordinates."""
from dataclasses import dataclass
import math

import numpy as np


@dataclass
class SecondaryLayout:
    coordinates: np.ndarray
    contexts: list[str]


def pair_table(length, pairs):
    """Reject repeated partners and crossings before constructing a planar tree."""
    partners = {}
    for i, j in pairs:
        if not 0 <= i < j < length or i in partners or j in partners:
            raise ValueError("Appariements incompatibles avec un dessin secondaire.")
        partners[i], partners[j] = j, i
    stack = []
    for i in range(length):
        if partners.get(i, -1) > i:
            stack.append(partners[i])
        elif i in partners:
            if not stack or stack.pop() != i:
                raise ValueError("Le dessin secondaire ne prend pas en charge les pseudonœuds.")
    return partners


def layout_secondary(length, pairs):
    """Draw stacked pairs as ladders and loop faces as regular polygons.

    Child helices grow radially out of their enclosing loop. Exterior components
    are separated by their bounding boxes. All distances are drawing units.
    An explicit work stack handles long helices without Python recursion.
    """
    if length < 1:
        raise ValueError("Une séquence non vide est requise.")
    partners = pair_table(length, pairs)
    coords = np.zeros((length, 2), dtype=float)
    contexts = ["Région non appariée"] * length
    spacing = 1.4
    cursor, index = 0., 0
    while index < length:
        end = partners.get(index, -1)
        if end < index:
            coords[index] = [cursor, 0.]
            cursor += spacing
            index += 1
            continue
        coords[index], coords[end] = [-spacing / 2, 0.], [spacing / 2, 0.]
        pending = [(index, end, np.array([0., 1.]))]
        while pending:
            i, j, outward = pending.pop()
            contexts[i] = contexts[j] = "Tige appariée"
            if partners.get(i + 1) == j - 1 and i + 1 < j - 1:
                coords[i + 1] = coords[i] + outward * spacing
                coords[j - 1] = coords[j] + outward * spacing
                pending.append((i + 1, j - 1, outward))
                continue
            vertices, children, unpaired = [i], [], []
            k = i + 1
            while k < j:
                other = partners.get(k, -1)
                if other > k:
                    children.append((k, other))
                    vertices.extend([k, other])
                    k = other + 1
                else:
                    vertices.append(k)
                    unpaired.append(k)
                    k += 1
            vertices.append(j)
            if len(vertices) <= 2:
                continue
            if not children:
                context = "Boucle terminale"
            elif len(children) > 1:
                context = "Jonction multibranche"
            else:
                a, b = children[0]
                context = "Boucle interne" if a > i + 1 and b < j - 1 else "Renflement"
            for k in unpaired:
                contexts[k] = context
            tangent = (coords[j] - coords[i]) / spacing
            angle = math.pi / len(vertices)
            radius = spacing / (2 * math.sin(angle))
            center = (coords[i] + coords[j]) / 2 + outward * spacing / (2 * math.tan(angle))
            for number, k in enumerate(vertices[1:-1], 1):
                theta = -math.pi / 2 - angle - number * 2 * angle
                coords[k] = center + radius * (math.cos(theta) * tangent + math.sin(theta) * outward)
            for a, b in children:
                direction = (coords[a] + coords[b]) / 2 - center
                direction /= np.linalg.norm(direction)
                pending.append((a, b, direction))
        component = coords[index:end + 1]
        component[:, 0] += cursor - component[:, 0].min()
        cursor = component[:, 0].max() + spacing * 2
        index = end + 1
    return SecondaryLayout(coords, contexts)


def comparison_layouts(reference, mutant):
    """Share coordinates only when both pair sets fit the same nested scaffold."""
    wt = layout_secondary(len(reference["sequence"]), reference["base_pairs"])
    if mutant is None:
        return wt, None, "Référence seule"
    mut = layout_secondary(len(mutant["sequence"]), mutant["base_pairs"])
    if len(reference["sequence"]) != len(mutant["sequence"]):
        return wt, mut, "Longueurs différentes · sans alignement"
    union = sorted(set(map(tuple, reference["base_pairs"])) | set(map(tuple, mutant["base_pairs"])))
    try:
        shared = layout_secondary(len(reference["sequence"]), union).coordinates
    except ValueError:
        # Changed partners can be incompatible: each panel must show its own topology.
        return wt, mut, "Dispositions propres à chaque structure"
    wt.coordinates = shared
    mut.coordinates = shared.copy()
    return wt, mut, "Positions fixées pour comparer WT et mutant"
