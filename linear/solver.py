import copy
from collections import deque
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--infile', help='File path to read matrix from')
parser.add_argument('--outfile', help='File path to write matrix variable values to')


class Matrix:
    def __init__(self, matrix):
        self.matrix = []
        for eq in matrix:
            self.matrix.append(list(map(lambda x: complex(x), eq)))

    @classmethod
    def from_file(cls, path):
        with open(path) as file:
            string_queue = deque((str_matrix := file.read()).split())
            print(str_matrix)
            if len(string_queue) < 4:
                raise ValueError('At least number of rows, number of columns '
                                 ', one variable and value for it must be in the file')
            vars_num = int(string_queue.popleft())
            eq_num = int(string_queue.popleft())
            if eq_num * (vars_num + 1) != len(string_queue):
                raise ValueError('Wrong number of elements')

            matrix = []
            for i in range(eq_num):
                matrix.append([])
                for _ in range(vars_num + 1):
                    matrix[i].append(string_queue.popleft())
            return cls(matrix)

    def solve(self):

        variables = []
        matrix = copy.deepcopy(self.matrix)
        col_permutations = deque()
        print('__________________________________________________')

        try:
            Matrix._check_zero_equations(matrix)
            i, j = 0, 0
            while i < len(matrix) and j < len(matrix[0]) - 1:
                if matrix[i][j] == 0:
                    Matrix._swap_zero_entry(matrix, i, j, col_permutations)

                if (coef := matrix[i][j]) != 1:
                    matrix[i] = list(map(lambda x: x / coef, matrix[i]))
                    print(f'R{i + 1} / {str(coef).strip("()")} -> R{i + 1}')

                Matrix._eliminate(matrix, i)
                Matrix._check_zero_equations(matrix)
                i, j = i + 1, j + 1

            i, j = i - 1, j - 1
            while i >= 0 and j >= 0:
                Matrix._eliminate(matrix, i, reverse=True)
                Matrix._check_zero_equations(matrix)
                i, j = i - 1, j - 1

            for _ in range(len(col_permutations)):
                Matrix._swap_cols(matrix, *reversed(col_permutations.pop()))

            for i, equation in enumerate(matrix):
                if len(list(filter(lambda x: x != 0, equation))) > 2:
                    raise Matrix.InfiniteSolutionsException
                variables.append(Matrix._get_formatted_complex(equation[-1]))

        except (Matrix.NoSolutionException, Matrix.InfiniteSolutionsException) as err:
            variables = [err]
        return variables

    class NoSolutionException(Exception):
        def __str__(self): return 'No solutions'

    class InfiniteSolutionsException(Exception):
        def __str__(self): return 'Infinitely many solutions'

    @staticmethod
    def _eliminate(matrix, entry_, reverse=False):
        start = entry_ + 1 if not reverse else entry_ - 1
        end = len(matrix) if not reverse else -1
        step = 1 if not reverse else -1

        for other in range(start, end, step):
            if (coef_ := -matrix[other][entry_]) == 0:
                continue
            matrix[other] = list(map(lambda x, y: x + y * coef_, matrix[other], matrix[entry_]))
            print(f'{str(coef_).strip("()")} * R{entry_ + 1} + R{other + 1} -> R{other + 1}')

    @staticmethod
    def _check_zero_equations(matrix):
        for eq_ in matrix:
            vars_ = set(eq_[:-1])
            if len(vars_) == 1 and vars_.pop() == 0:
                if eq_[-1] != 0:
                    raise Matrix.NoSolutionException
                else:
                    matrix.remove(eq_)

    @staticmethod
    def _swap_zero_entry(matrix, i_, j_, perm_stack):
        for col in range(j_, len(matrix[0]) - 1):
            for row in range(i_, len(matrix)):
                if matrix[row][col] != 0:
                    if row != i_:
                        print(f'R{i_ + 1} <-> R{row + 1}')
                        Matrix._swap_rows(matrix, i_, row)
                    if col != j_:
                        print(f'C{j_ + 1} <-> C{col + 1}')
                        Matrix._swap_cols(matrix, j_, col)
                        perm_stack.append((j_, col))
                    return

    @staticmethod
    def _swap_rows(matrix, row_a, row_b):
        matrix[row_a], matrix[row_b] = matrix[row_b], matrix[row_a]

    @staticmethod
    def _swap_cols(matrix, col_a, col_b):
        for row_ind in range(len(matrix)):
            matrix[row_ind][col_a], matrix[row_ind][col_b] = matrix[row_ind][col_b], matrix[row_ind][col_a]

    @staticmethod
    def _get_formatted_complex(num: complex):
        return num if num.imag != 0 else num.real


if __name__ == '__main__':
    args = parser.parse_args()
    with open(args.outfile, 'w') as out_file:
        input_matrix = Matrix.from_file(args.infile)
        result = list(map(lambda x: str(x).strip('()'), input_matrix.solve()))
        out_file.write('\n'.join(result))
        print('__________________________________________________', *result, sep='\n')
