"""A collection of miscellaneous functions used for Cluster Expansion."""
from itertools import permutations, combinations, product
import numpy as np
import collections


def index_by_position(atoms):
    """Set atomic indices by its position."""
    # add zero to avoid negative zeros
    tags = atoms.get_positions().round(decimals=6) + 0
    tags = tags.tolist()
    deco = sorted([(tag, i) for i, tag in enumerate(tags)])
    indices = [i for tag, i in deco]
    return indices


def sort_by_position(atoms):
    """Sort atoms by its position."""
    # Return a new Atoms object with sorted atomic order.
    # The default is to order according to chemical symbols,
    # but if *tags* is not None, it will be used instead.
    # A stable sorting algorithm is used.
    indices = index_by_position(atoms)
    return atoms[indices]


def wrap_and_sort_by_position(atoms):
    """Wrap and sort atoms by their positions."""
    atoms.wrap()
    atoms = sort_by_position(atoms)
    return atoms


def create_cluster(atoms, indices):
    """Create a cluster centered in the unit cell."""
    cluster = atoms[list(indices)]
    cell = cluster.get_cell()
    center = 0.5 * (cell[0, :] + cell[1, :] + cell[2, :])
    min_max_dist = 1E10
    minimal_cluster = None
    for ref_indx in range(len(indices)):
        cluster_cpy = cluster.copy()
        sub_indx = [i for i in range(len(indices)) if i != ref_indx]
        mic_dists = cluster_cpy.get_distances(
            ref_indx, sub_indx, mic=True, vector=True)
        com = cluster_cpy[ref_indx].position + \
            np.sum(mic_dists, axis=0) / len(indices)
        cluster_cpy.translate(center - com)
        cluster_cpy.wrap()
        pos = cluster_cpy.get_positions()
        lengths = []

        for comb in combinations(range(len(indices)), r=2):
            dist = pos[comb[0], :] - pos[comb[1], :]
            length = np.sqrt(np.sum(dist**2))
            lengths.append(length)
        max_dist = np.max(lengths)
        if max_dist < min_max_dist:
            min_max_dist = max_dist
            minimal_cluster = cluster_cpy.copy()
    return minimal_cluster


def shift(array):
    ref = array[-1]
    array[1:] = array[:-1]
    array[0] = ref
    return array


def distances_and_angles(atoms, ref_indx, float_obj_angle):
    """Get sorted internal angles of a."""
    indices = [a.index for a in atoms if a.index != ref_indx]
    if len(atoms) < 2:
        raise ValueError("distances and angles cannot be called for"
                         "{} body clusters".format(len(atoms)))
    if len(atoms) == 2:
        dist = atoms.get_distance(ref_indx, indices[0])
        classifier = distance_string(atoms.info["distances"], dist)
        return [classifier]

    angles = []
    dists = []

    for comb in combinations(indices, r=2):
        angle = atoms.get_angle(comb[0], ref_indx, comb[1])
        ang_classifier = float_obj_angle.get(angle)
        angles.append(ang_classifier)

    dists = atoms.get_distances(ref_indx, indices, mic=True)
    dists = sorted(dists.tolist(), reverse=True)
    dists = [distance_string(atoms.info["distances"], d) for d in dists]
    return dists + sorted(angles, reverse=True)


def get_cluster_descriptor(cluster, float_obj_angle):
    """Create a unique descriptor for each cluster."""
    dist_ang_tuples = []
    for ref_indx in range(len(cluster)):
        dist_ang_list = distances_and_angles(cluster, ref_indx,
                                             float_obj_angle)
        dist_ang_tuples.append(dist_ang_list)
    return dist_ang_tuples


def sort_by_internal_distances(atoms, indices, float_obj_ang):
    """Sort the indices according to the distance to the other elements."""
    if len(indices) <= 1:
        return list(range(len(indices))), "point"

    cluster = create_cluster(atoms, indices)
    if len(indices) == 2:
        dist_ang = get_cluster_descriptor(cluster, float_obj_ang)
        order = list(range(len(indices)))
        eq_sites = [(0, 1)]
        descr = "{}_0".format(dist_ang[0][0])
        return order, eq_sites, descr

    dist_ang = get_cluster_descriptor(cluster, float_obj_ang)
    sort_order = [ind for _, ind in sorted(zip(dist_ang, range(len(indices))))]
    dist_ang.sort()
    equivalent_sites = [[i] for i in range(len(indices))]
    site_types = [i for i in range(len(indices))]
    for i in range(len(sort_order)):
        for j in range(i + 1, len(sort_order)):
            if dist_ang[i] == dist_ang[j]:
                if site_types[j] > i:
                    # This site has not been assigned to another category yet
                    site_types[j] = i
                st = site_types[j]
                if j not in equivalent_sites[st]:
                    equivalent_sites[st].append(j)

    # Remove empty lists from equivalent_sites
    equivalent_sites = [entry for entry in equivalent_sites if len(entry) > 1]

    # Create a string descriptor of the clusters
    dist_ang_strings = []
    for item in dist_ang:
        strings = [str(x) for x in item]
        dist_ang_strings.append("_".join(strings))
    string_description = "-".join(dist_ang_strings)
    return sort_order, equivalent_sites, string_description


def ndarray2list(data):
    """
    Convert nested lists of a combination of lists and numpy arrays
    to list of lists
    """
    if not isinstance(data, list) and not isinstance(data, np.ndarray):
        return data

    data = list(data)
    for i in range(len(data)):
        data[i] = ndarray2list(data[i])
    return list(data)


def dec_string(deco, equiv_sites):
    """Create the decoration string based on equiv sites."""
    equiv_dec = sorted(equivalent_deco(deco, equiv_sites))
    return ''.join(str(i) for i in equiv_dec[0])


def equivalent_deco(deco, equiv_sites):
    """Generate equivalent decoration numbers based on equivalent sites."""
    if not equiv_sites:
        return [deco]

    perm = []
    for equiv in equiv_sites:
        perm.append(list(permutations(equiv)))

    equiv_deco = []
    for comb in product(*perm):
        order = []
        for item in comb:
            order += list(item)

        orig_order = list(range(len(deco)))
        for i, srt_indx in enumerate(sorted(order)):
            orig_order[srt_indx] = order[i]
        equiv_deco.append([deco[indx] for indx in orig_order])

    unique_deco = []
    for eq_dec in equiv_deco:
        if eq_dec not in unique_deco:
            unique_deco.append(eq_dec)
    return unique_deco


def flatten(x):
    """Flatten list."""
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]


def get_unique_name(size, max_dia, fam_id):
    """Get unique cluster names."""
    name = "c{}_{}_{}".format(size, max_dia, fam_id)
    return name


def nested_array2list(array):
    """Convert a nested array/tuple to a nested list."""
    if isinstance(array, np.ndarray):
        array = array.tolist()
    else:
        array = list(array)
    try:
        for i in range(len(array)):
            if isinstance(array[i], np.ndarray):
                array[i] = array[i].tolist()
            else:
                array[i] = list(array[i])
    except TypeError:
        pass
    return array


def update_db(uid_initial=None, final_struct=None, db_name=None,
              custom_kvp_init={}, custom_kvp_final={}):
    """Update the database.

    Arguments:
    =========
    uid_initial: int
        entry ID of the initial structure in the database

    final_struct: Atoms
        Atoms object with the final structure with a physical
        quantity that needs to be modeled (e.g., DFT energy)

    db_name: str
        Database name

    custom_kvp_init: dict (optional)
        If desired, one can pass additional key-value-pairs for the
        entry containing the initial structure

    custom_kvp_final: dict (optional)
        If desired, one can pass additional key-value-pairs for the
        entry containing the final structure
    """
    from ase.db import connect
    db = connect(db_name)

    init_row = db.get(id=uid_initial)

    # Check if a final structure already exits
    name = init_row.name
    select_cond = [('name', '=', name), ('struct_type', '=', 'final')]
    exist = sum(1 for row in db.select(select_cond))
    if exist >= 1:
        print("A structure with 'name'={} and 'struct_type'=final "
              "already exits in DB".format(name))
        return

    # Write the final structure to database
    kvp_final = {'struct_type': 'final', 'name': name}
    kvp_final.update(custom_kvp_final)
    uid_final = db.write(final_struct, key_value_pairs=kvp_final)

    kvp_update_init = {
        'converged': True,
        'started': '',
        'queued': '',
    }
    if kvp_final['struct_type'] == 'final':
        kvp_update_init['final_struct_id'] = uid_final
    kvp_update_init.update(custom_kvp_init)

    # Update info for the initial structure
    db.update(uid_initial, **kvp_update_init)


def exclude_information_entries():
    """Return selection condition to exlcude all entries in the database that
       only contain information about the clusters.
    """
    return [('name', '!=', 'unit_cell'),
            ('name', '!=', 'template'),
            ('name', '!=', 'float_classification')]


def get_all_internal_distances(atoms, max_dist):
    """Obtain all internal distances of the passed atoms object and return a
       Numpy array containing all the distances sorted in an ascending order.
    """
    from scipy.spatial import cKDTree as KDTree

    tree = KDTree(atoms.get_positions())
    distances = []
    for atom in atoms:
        indices = tree.query_ball_point(atom.position, max_dist)
        dists = atoms.get_distances(atom.index, indices)
        for d in dists:
            if np.any(np.abs(np.array(distances) - d) < 1E-6):
                continue
            distances.append(d)
    distances = sorted(distances)

    # Make sure that the first element is 0
    assert distances[0] < 1E-6
    return np.array(distances[1:])


def distance_string(distance_array, distance):
    """Provide a name of the passed distance in terms of the nearest neighbor
       based on the internal distances in an array."""
    indx = np.argmin(np.abs(distance_array - distance))
    assert abs(distance_array[indx] - distance) < 1E-6

    if indx < 9:
        return "0{}nn".format(indx+1)
    return "{}nn".format(indx+1)


def reconfigure(setting, select_cond=None):
    from ase.clease import CorrFunction
    setting.reconfigure_settings()
    CorrFunction(setting).reconfigure_db_entries(select_cond)


def split_dataset(X, y, nsplits=10):
    """Split the dataset such that it can be used for k-fold
        cross validation."""
    from random import shuffle
    indices = list(range(len(y)))
    shuffle(indices)
    partitions = []
    num_validation = int(len(y)/nsplits)

    if num_validation < 1:
        num_validation = 1
    for i in range(nsplits):
        start = i*num_validation
        end = (i+1)*num_validation
        if i == nsplits-1:
            indx = indices[start:]
        else:
            indx = indices[start:end]
        mask = np.zeros(len(y), dtype=np.uint8)
        mask[indx] = 1
        data = {
            "train_X": X[mask == 0, :],
            "train_y": y[mask == 0],
            "validate_X": X[mask == 1, :],
            "validate_y": y[mask == 1]
        }
        partitions.append(data)
    return partitions
