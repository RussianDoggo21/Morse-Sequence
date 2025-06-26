from simplextree import SimplexTree as SimplexTreeCpp
from morse_sequence._core import MorseSequence as MScpp

import numpy as np
from typing import Collection, Iterable

class MorseSequence(MScpp):
  def __init__(self, st: SimplexTreeCpp) -> None:
    # Call the C++ constructor (via pybind11)
    super().__init__(st)

    # Locally retain the SimplexTree (optional but useful)
    self.st = st


	  