from simplextree import SimplexTree
import morse_sequence._core as _core

def check_simplextree_conflict():

    print("=== SimplexTree from simplextree ===")
    print(SimplexTree, SimplexTree.__module__)

    # Créons un SimplexTree Python
    st = SimplexTree([[1, 2, 3]])

    # Tentative de création d'un MorseSequence
    try:
        ms = _core.MorseSequence(st)
        print("\n✅ Pas de conflit : le SimplexTree de simplextree est accepté par MorseSequence.")
    except TypeError as e:
        print("\n❌ Conflit détecté : le SimplexTree de simplextree n'est pas reconnu par MorseSequence.")
        print("Message d'erreur :", e)

    # Vérifions aussi le MRO (Method Resolution Order) de MorseSequence
    print("\n=== MRO (Method Resolution Order) of MorseSequence ===")
    for t in _core.MorseSequence.__mro__:
        print(" -", t, getattr(t, "__module__", "?"))

    # Affichons la classe SimplexTree
    print("\n=== SimplexTree class object ===")
    print(SimplexTree)


# Lancer directement le check
#check_simplextree_conflict()


from simplextree._simplextree import SimplexTree
import morse_sequence._core as core

st = SimplexTree([[1,2,3]])
ms = core.MorseSequence(st)  # Ça fonctionne car c'est le bon type


"""
st = SimplexTree([[1, 2, 3]])
print(type(st))                    # devrait être <class 'simplextree.SimplexTree.SimplexTree'>
print(st.__class__.__module__)     # devrait être 'simplextree.SimplexTree'
print(_core.MorseSequence.__init__.__doc__)  # montre le type attendu
"""

"""
import simplextree._simplextree as st1
import morse_sequence._core as core

print(st1.SimplexTree is core.__dict__['SimplexTree'])
"""

"""
from simplextree import _simplextree
import morse_sequence._core as core

st1 = _simplextree.SimplexTree([[1,2,3]])
st2 = _simplextree.SimplexTree([[1,2,3]])

print(type(st1))  # type C++ pybind11
print(type(st2))

print(isinstance(st1, core.MorseSequence.__init__.__annotations__['st']))
"""

# To run the file from the root : python3 -m tests.python_tests.test_version