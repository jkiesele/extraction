# -*- coding: utf-8 -*-
"""
A luigi task to fill histograms
"""
import luigi

from workflow.BaseTask import BaseTask
from workflow.HistoTask import AllHistoTasks

from config.histograms import histograms
from config.regions import regions
from config.systematics import systematics


class BranchPlotTask(BaseTask):
    year = luigi.IntParameter()
    region = luigi.Parameter()
    systematic = luigi.Parameter()
    histogram = luigi.Parameter()

    def log(self):
        return f'./logs/branch_plot/{self.year}/{self.region}/{self.systematic}/{self.histogram}.log'

    def requires(self):
        return [AllHistoTasks(year=self.year)]

    def output(self):
        return [luigi.LocalTarget(self.log())]

    def run(self):
        self.save_execute(command=f'python python/plot_branch.py --year {self.year} --region {self.region} \
                                    --systematic {self.systematic} --histo {self.histogram}', log=self.log())



class AllBranchPlotTasks(luigi.WrapperTask):
    year = luigi.IntParameter()

    def requires(self):
        for histogram in histograms.keys():
            for region in regions.keys():
                for systematic in systematics.keys():
                    yield BranchPlotTask(year=self.year, region=region, systematic=systematic, histogram=histogram)