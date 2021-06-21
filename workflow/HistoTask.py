# -*- coding: utf-8 -*-
"""
A luigi task to fill histograms
"""
import luigi

from workflow.BaseTask import BaseTask

from config.general import histopath
from config.samples import samples
from config.regions import regions
from config.systematics import systematics


class HistoTask(BaseTask):
    year = luigi.IntParameter()
    sample = luigi.Parameter()
    region = luigi.Parameter()
    systematic = luigi.Parameter()

    def log(self):
        return f'./logs/fill_histos/{self.year}/{self.region}/{self.systematic}/{self.sample}.log'

    def output(self):
        return [luigi.LocalTarget(self.log()), luigi.LocalTarget(histopath(isMC=samples[self.sample]['MC'],
                                                                           year=self.year,
                                                                           filename=self.sample,
                                                                           region=self.region,
                                                                           systematic=self.systematic))]

    def run(self):
        self.save_execute(command=f'python python/fill_histos.py --year {self.year} --sample {self.sample} \
                                    --region {self.region} --systematic {self.systematic}', log=self.log())



class AllHistoTasks(luigi.WrapperTask):
    year = luigi.IntParameter()

    def requires(self):
        for sample in samples.keys():
            for region in regions.keys():
                for systematic in systematics.keys():
                    yield HistoTask(year=self.year, sample=sample, region=region, systematic=systematic)