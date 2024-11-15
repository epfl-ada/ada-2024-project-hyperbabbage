from abc import ABC, abstractmethod
import re
from typing import List
import pandas as pd


"""
This class represents the data needed to convert a column in the raw dataframe, to a column in the new dataframe.
"""
class ColumnClean:
    def __init__(self, old_column_name: str, new_column_name: str, clean = lambda x: x) -> None:
        self.old_column_name = old_column_name
        self.new_column_name = new_column_name
        self.clean = clean
    def run(self, old_df, new_df):
        new_df[self.new_column_name] = old_df[self.old_column_name].apply(lambda x: self.clean(x))


class Preprocessing:
    def __init__(self, strategies: List[ColumnClean]) -> None:
        self.strategies = strategies
    
    def transform(self, old_df: pd.DataFrame) -> pd.DataFrame:
        new_df = pd.DataFrame(index=old_df.index)
        for stragegy in self.strategies:
            stragegy.run(old_df, new_df)
        return new_df.convert_dtypes()
    
    def get_all_old_columns(self):
        return [strategy.old_column_name for strategy in self.strategies]

    def get_all_new_columns(self):
        return [strategy.new_column_name for strategy in self.strategies]
