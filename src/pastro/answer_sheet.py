from pathlib import Path
from typing import Any, Callable
import yaml
import inspect
import re


class AnswerSheet:
    """Base class for answer sheet generators

    :param answer_sheet: Path to answer sheet. Default is 'answer-sheet.yaml'.
    :param readme: Path to README.md as string. Default is 'README.md'.
    """

    def __init__(self, answer_sheet: str = "answer-sheet.yaml",
                 readme: str = "README.md") -> None:

        self.answer_sheet = Path(answer_sheet)
        self.__get_expected_answers(readme)

        return None

    def __get_expected_answers(self, readme: str) -> None:
        """Retrieve expected answers for statement

        Fetch answer keys from the assignment's README and save them as keys
        of the `self.answers` dictionary. Values are initialized to `None`.

        :param readme: Path to README.md file.
        """

        # TODO: Add test to ensure that values of questions add up to 10

        # Define REGEX for answer keys
        p = re.compile(
            r"""[sl][fibs]_[RMAC][ED][GTX]_[^_]+_[^_]+_[\d.]+""")

        # Initialize key list
        self.answers = {}

        # Read statement and retrieve keys
        with open(readme) as source:
            for line in source:
                keys = p.findall(line)
                for item in keys:
                    self.answers[item] = None

        return None

    def add_answer(self, answer: str, value: Any) -> None:
        """DEPRECATED: USE self.add INSTEAD"""

        if answer in self.answers.keys():
            self.answers[answer] = value
        else:
            raise KeyError(f"Answer '{answer}' not in expected answers.")

        return None

    def add(self, answer: str, value: Any) -> None:
        """Setter for answers dictionary.

        Adds answer to dictionary if the key is already present in
        self.answers. Raises KeyError otherwise.

        :param answer: Answer key
        :param value: Answer value
        """

        # TODO: Add type checker according to code in answer key

        if answer in self.answers.keys():
            self.answers[answer] = value
        else:
            raise KeyError(f"Answer '{answer}' not in expected answers.")

        return None

    def __call__(self) -> None:
        """Populate answer sheet with data from answers dictionary."""

        # Give warning if any of the answers is empty
        for key, val in self.answers.items():
            if val is None:
                print(f"Question {key} is not answered!")

        with open(self.answer_sheet, "w") as file:
            yaml.dump(self.answers, file)

        return None


def report_code(functions: Callable | list[Callable]) -> str:
    """Create string with filename and source code for given functions.

    :param functions: List of functions to be reported.
    :return: Filename and source code.
    """

    if isinstance(functions, Callable):
        functions = [functions]

    # Get path to file relative to CWD
    filenames = [Path(inspect.getfile(function)).relative_to(Path.cwd())
                 for function in functions]

    # Get source code
    sources = [inspect.getsource(function) for function in functions]

    return "\n\n".join([f"{filename}:\n{source}"
                        for filename, source in zip(filenames, sources)])


def report_function(function: Callable) -> str:
    """Create string with filename and source code for given function.

    NOTE: DEPRECATED. Use `report_code` instead.

    :param function: Function to be reported.
    :return: Filename and source code.
    """

    # Get path to file relative to CWD
    filename = Path(inspect.getfile(function)).relative_to(Path.cwd())

    # Get source code
    source = inspect.getsource(function)

    return f"{filename}:\n{source}"
