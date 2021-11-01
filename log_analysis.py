from abc import ABCMeta, abstractmethod
import re
import numpy as np
import pandas as pd


class LogMain(metaclass=ABCMeta):
    """
    Abstract object
    All classes that will inherit from this class will need to include the following methods
    """

    @abstractmethod
    def read_log(self):
        """
        Method to read log and provide as a file type which is most according depending on the .log file
        """
        raise NotImplementedError()

    @abstractmethod
    def check_log(self, *args, **kwargs):
        """
        Perform checks on the log. Ideally each class will have other methods which will help this method out
        """
        raise NotImplementedError()


class Bwa(LogMain):
    """
    This class will check the bwa log
    """
    # TODO: Can we somehow include the numbers that the log provides us for anything? (e.g. percentiles, mean, std, ...)
    # TODO: What is the threshold for "enough pairs"?
    # TODO: Check if the output expected actually exists

    def __init__(self, path, sample):
        self.log_file = None
        self.path = path
        self.sample = sample
        self.dict_ = None
        self.paired = self.single_paired()
        self.read_log()

    def single_paired(self, table_path='data/fastq.csv'):
        """
        Check whether the sample is paired (R1 and R2) or single (only R1). This method is relevant
        for cases in which two different log files are generated
        :param table_path: Path in which we can find the fastq.csv (file containing this information)
        :return: Boolean value which will be stored as part of the class variables
        """

        df = pd.read_csv(table_path)

        bool_ = len(df[df.Sample == self.sample]) == 2

        return bool_

    def read_log(self):
        """
        Method to store the log file as part of the class variables
        """
        with open(self.path + self.sample + '.log') as f:
            self.log_file = f.readlines()

    def check_log(self):
        """
        This method will run all the methods implemented for this class

        - ``check_lines()``
        - ``check_num_sequence()``
        - ``check_consistency()``
        - ``check_start_statement()``
        - ``check_correct_sample()``
        - ``check_candidate_pairs()``
        - ``check_not_enough_pairs(threshold=100)``
        - ``check_enough_pairs(threshold=100)``
        """
        self.check_lines()
        self.check_num_sequence()
        self.check_consistency()
        self.check_start_statement()
        self.check_correct_sample()
        self.check_candidate_pairs()
        self.check_not_enough_pairs(threshold=100)
        self.check_enough_pairs(threshold=100)

    def check_lines(self):
        """
        Check the number of lines of the log:

        - In case of paired samples (R1 + R2) the log should contain 16 lines
        - In case of single samples (R1) the log should contain 6 lines
        """
        if self.paired:
            if len(self.log_file) != 16:
                raise Exception('check_lines: ' + self.sample + ' has the wrong number of log lines')
        else:
            if len(self.log_file) != 6:
                raise Exception('check_lines: ' + self.sample + ' has the wrong number of log lines')

    def check_num_sequence(self):
        """
        Make sure we are reading a positive number of sequences

        - ``[M::process] read 400000 sequences (40000000 bp)...``
        - 400000 must be a positive number
        """
        seq = re.findall(r'\d+', self.log_file[1])
        if any(int(i) < 0 for i in seq):
            raise Exception('check_num_sequence: ' + self.sample + ' did not read a positive number of sequences')

    def check_consistency(self):
        """
        Check consistency in terms of single-end and paired-end sequences. This check is only done for paired samples

        - ``[M::process] read 400000 sequences (40000000 bp)...``
        - ``[M::process] 0 single-end sequences; 400000 paired-end sequences``
        - 40000 (in the first line) has to be equal to 0 + 400000 (second line)
        """
        if self.paired:
            seq1 = int(re.findall(r'\d+', self.log_file[1])[0])
            seq2 = list(map(int, re.findall(r'\d+', self.log_file[2])))
            if seq1 != sum(seq2):
                raise Exception('check_consistency: ' + self.sample + ' has inconsistency in terms of paired-end and '
                                                                      'single-end sequences')
        else:
            seq1 = int(re.findall(r'\d+', self.log_file[1])[0])
            seq2 = int(re.findall(r'\d+', self.log_file[2])[0])
            if seq1 != seq2:
                raise Exception('check_consistency: ' + self.sample + ' has inconsistency in terms of read sequences '
                                                                      'and processed sequences')

    def check_start_statement(self):
        """
        Make sure that the initial log statement is the expected one

        - ``[M::bwa_idx_load_from_disk] read 3171 ALT contigs``
        """
        if self.log_file[0][:-1] != '[M::bwa_idx_load_from_disk] read 3171 ALT contigs':
            raise Exception('check_start_statement: ' + self.sample + ' does not have the correct log starting '
                                                                      'statement')

    def check_finish_statement(self):
        """
        Make sure that the final log statement is the expected one. We check that the text is the expected one and also
        that the runtime is a positive value

        - ``[main] Real time: 23.460 sec; CPU: 133.799 sec``
        - Check for [main] Real time:
        - Check real time and CPU times are both positive
        """
        text = self.log_file[-1][:17]
        nums = list(map(int, re.findall(r"[-+]?\d*\.\d+|\d+", self.log_file[-1])))
        if any(i <= 0 for i in nums) | text != '[main] Real time:':
            raise Exception('check_finish_statement: ' + self.sample + ' does not have the final statement we expected')

    def check_correct_sample(self):
        """
        Make sure that the file that is being processed is the actual sample

        - ``[main] CMD: bwa mem -p -t 8 -R @RG\tID:HSRR062625\tLB:HSRR062625\tSM:HSRR062625\tPU:unknown\ ...``
        - If we are processing sample HSRR062625 we should only have this value in this string
        """
        if len(re.findall(self.sample, self.log_file[-2])) == 0:
            raise Exception('check_correct_sample: ' + self.sample + ' should be processed however another sample has '
                                                                     'been processed instead')

    def check_candidate_pairs(self):
        """
        Check the candidate pairs are being evaluated. This check is only done for paired samples

        - ``[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR):``
        - We make sure that this string is present in the log file
        """
        if self.paired:
            if '[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR):' != self.log_file[3][:62]:
                raise Exception(
                    'check_candidate_pairs: ' + self.sample + ' does not have the expected candidate pairs output')

    def check_not_enough_pairs(self, threshold=100):
        """
        Whenever there are not enough pairs the command skips those pairs. This check is only done for paired samples

        - ``[M::mem_pestat] skip orientation FF as there are not enough pairs``
        - We make sure that the orientation FF actually does not have enough pairs

        :param threshold: Threshold used to define whether we actually have enough pairs
        """
        if self.paired:
            FF, FR, RF, RR = list(map(int, self.log_file[3][64:-2].strip().split(',')))
            self.dict_ = {'FF':FF,
                     'FR': FR,
                     'RF': RF,
                     'RR': RR}
            for i in self.log_file[4:-4]:
                if i[:32] == '[M::mem_pestat] skip orientation':
                    if self.dict_[i[33:35]] > threshold:
                        raise Exception('check_not_enough_pairs ' + self.sample + ' has skipped an orientation due to low '
                                                                              'number of pairs where this is not the '
                                                                              'case')

    def check_enough_pairs(self, threshold=100):
        """
        When we have enough pairs, there are several things that need to be checked. This check is only done for paired
        samples

        - ``[M::mem_pestat] analyzing insert size distribution for orientation FR...``
        - This log line is only shown for orientations with sufficient pairs (FR in this case). Once there are enough pairs 4 steps are done on them. We check that these steps are really done

        :param threshold: Threshold used to define whether we actually have enough pairs
        """
        if self.paired:
            for num, i in enumerate(self.log_file[4:-4]):
                if i[-4:-1] == '...':
                    if self.dict_[self.log_file[4:-4][1][-6:-4]] > threshold:
                        if ((self.log_file[4 + num + 1][:40] != '[M::mem_pestat] (25, 50, 75) percentile:') |
                                (self.log_file[4 + num + 2][
                                 :71] != '[M::mem_pestat] low and high boundaries for computing mean and std.dev:') |
                                (self.log_file[4 + num + 3][:33] != '[M::mem_pestat] mean and std.dev:') |
                                (self.log_file[4 + num + 4][
                                 :57] != '[M::mem_pestat] low and high boundaries for proper pairs:')):
                            raise Exception('check_enough_pairs ' + self.sample + ' not all the steps of BWA have been '
                                                                                  'executed')
                    else:
                        raise Exception('check_enough_pairs ' + self.sample + ' not enough pairs, however code has '
                                                                              'considered to have enough pairs anyway')


class Fastqc(LogMain):
    """
    This class will check the fastqc log
    """

    def __init__(self, path, sample):
        self.log_file_1 = None
        self.log_file_2 = None
        self.path = path
        self.sample = sample
        self.read_log()

    def read_log(self):
        """
        Method to store the log file as part of the class variables
        """
        try:
            with open(self.path + self.sample + '_R1_fastqc.log') as f:
                self.log_file_1 = f.readlines()
        except:
            self.log_file_1 = False

        try:
            with open(self.path + self.sample + '_R2_fastqc.log') as f:
                self.log_file_2 = f.readlines()
        except:
            self.log_file_2 = False

    def check_log(self):
        """
        This method will run all the methods implemented for this class

        - ``check_lines()``
        - ``check_start_end()``
        """
        self.check_lines()
        self.check_start_end()

    def check_lines(self):
        """
        Check correct number of lines in log

        - We expect 22 lines in the log
        """
        if self.log_file_1:
            if len(self.log_file_1) != 22:
                raise Exception('check_lines: ' + self.sample + '_R1 does not have the correct number of lines')

        if self.log_file_2:
            if len(self.log_file_2) != 22:
                raise Exception('check_lines: ' + self.sample + '_R2 does not have the correct number of lines')

    def check_start_end(self):
        """
        Check start and end statements are the expected ones

        - ``Started analysis`` is the start log line
        - ``Analysis complete`` is the end log line
        """
        if self.log_file_1:
            if (self.log_file_1[0][:16] != 'Started analysis') | (self.log_file_1[-1][:17] != 'Analysis complete'):
                raise Exception('check_lines: ' + self.sample + '_R1 does not seem to have been processed properly')

        if self.log_file_2:
            if (self.log_file_2[0][:16] != 'Started analysis') | (self.log_file_2[-1][:17] != 'Analysis complete'):
                raise Exception('check_lines: ' + self.sample + '_R2 does not seem to have been processed properly')


class SamSort(LogMain):
    """
    This class will check the samsort log
    """
    # TODO: Make sure check_unmated makes sense (I think the mate them later on so we should expect them to be 0 anyway
    #  but not sure)

    def __init__(self, path, sample):
        self.log_file = None
        self.path = path
        self.sample = sample
        self.paired = self.single_paired()
        self.read_log()
        self.dups = None

    def single_paired(self, table_path='data/fastq.csv'):
        """
        Check whether the sample is paired (R1 and R2) or single (only R1). This method is relevant
        for cases in which two different log files are generated
        :param table_path: Path in which we can find the fastq.csv (file containing this information)
        :return: Boolean value which will be stored as part of the class variables
        """

        df = pd.read_csv(table_path)

        bool_ = len(df[df.Sample == self.sample]) == 2

        return bool_

    def read_log(self):
        """
        Method to store the log file as part of the class variables
        """
        with open(self.path + self.sample + '_sort_nodup.sam.log') as f:
            self.log_file = f.readlines()

    def check_log(self):
        """
        This method will run all the methods implemented for this class

        - ``check_lines()``
        - ``check_start_statement()``
        - ``check_finish_statement()``
        - ``check_correct_sample()``
        - ``check_third_line()``
        - ``check_unmated()``
        - ``check_header()``
        - ``check_rows()``
        - ``check_table_sums()``
        - ``check_removals()``
        """
        self.check_lines()
        self.check_start_statement()
        self.check_finish_statement()
        self.check_correct_sample()
        self.check_third_line()
        self.check_unmated()
        self.check_header()
        self.check_rows()
        self.check_table_sums()
        self.check_removals()

    def check_lines(self):
        """
        Check the number of lines of the log:

        - In case of paired samples (R1 + R2) the log should contain 15 lines
        - In case of single samples (R1) the log should contain 14 lines
        """
        if self.paired:
            if len(self.log_file) != 15:
                raise Exception('check_lines: ' + self.sample + ' has the wrong number of log lines')
        else:
            if len(self.log_file) != 14:
                raise Exception('check_lines: ' + self.sample + ' has the wrong number of log lines')

    def check_start_statement(self):
        """
        Make sure that the initial log statement is the expected one

        - ``samblaster: Version 0.1.26``
        """
        if self.log_file[0][:-1] != 'samblaster: Version 0.1.26':
            raise Exception('check_start_statement: ' + self.sample + ' does not have the correct log starting '
                                                                      'statement')

    def check_finish_statement(self):
        """
        Make sure that the final log statement is the expected one.

        - ``[bam_sort_core] merging from files and in-memory blocks...``
        - We remove the numbers from the string we check since they are variable
        """
        if re.sub(" \d+", "", self.log_file[-1][:-1]) != '[bam_sort_core] merging from files and in-memory blocks...':
            raise Exception('check_finish_statement: ' + self.sample + ' does not have the final statement we expected')

    def check_correct_sample(self):
        """
        Make sure that the file that is being processed is the actual sample

        - ``samblaster: Opening OUTPUT_NEW/HSRR062650/bwa/HSRR062650.sam for read``
        - If we are processing sample HSRR062650 we should only have this value in this string
        """
        if len(re.findall(self.sample, self.log_file[1][:-1])) == 0:
            raise Exception('check_correct_sample: ' + self.sample + ' should be processed however another sample has '
                                                                     'been processed instead')

    def check_third_line(self):
        """
        The third line of these logs are all constant:

        - ``samblaster: Outputting to stdout``
        """
        if self.log_file[2][:-1] != 'samblaster: Outputting to stdout':
            raise Exception('check_third_line: ' + self.sample + ' does not have the expected output in line 3')

    def check_unmated(self):
        """
        Make sure that we have 0 mated pairs

        - ``samblaster: Found 0 of 200000 (0.000%) total read ids are marked paired yet are unmated.``
        - We check 0 and 0.000
        """
        if ((int(re.findall(r'\d+', self.log_file[4])[0]) != 0) |
            (float(re.findall(r"[-+]?\d*\.\d+|\d+", self.log_file[4])[-1]) != float(0))):
            raise Exception('check_unmated: ' + self.sample + ' has mated pairs when it should not have them yet')

    def check_header(self):
        """
        Check that the header of the table is correct
        """
        if self.log_file[6].replace(" ", "")[:-1] != 'samblaster:PairTypeType_ID_Count%Type/All_IDsDup_ID_Count%Dups/Type_ID_Count%Dups/All_Dups%Dups/All_IDs':
            raise Exception('check_header: ' + self.sample + ' does not have the expected table header')


    @staticmethod
    def _check_digit(x):
        """
        Aux function, not relevant
        """
        dig = False
        try:
            float(x)
            dig = True
        except:
            pass
        return dig

    def check_rows(self):
        """
        Check that the row names are all present

        - In the case of paired it should be ['Both Unmapped', 'Orphan/Singleton', 'Both', 'Total']
        - In the case of singles it should be ['Unmapped Orphan/Singleton', 'Mapped Orphan/Singleton', 'Total']
        """
        # Extract Data
        list_ = []
        for i in self.log_file[8:-3]:
            x = [j for j in i[12:-1].split() if not self._check_digit(j)]
            list_.append(' '.join(x))

        # Check values
        if self.paired:
            if list_ != ['Both Unmapped', 'Orphan/Singleton', 'Both Mapped', 'Total']:
                raise Exception('check_rows: ' + self.sample + ' does not have the expected row names')
        else:
            if list_ != ['Unmapped Orphan/Singleton', 'Mapped Orphan/Singleton', 'Total']:
                raise Exception('check_rows: ' + self.sample + ' does not have the expected row names')

    def check_table_sums(self):
        """
        Check that the table provided by the log file is correct:

        - We are checking that the values and the sums add up
        """
        # Extract data from table
        list_ = []
        for i in self.log_file[8:-3]:
            x = [float(j) for j in i[12:-1].split() if self._check_digit(j)]
            list_.append(x)

        self.dups = int(list_[-1][2])

        # Perform the summation accordingly
        sum_ = lambda x: np.isclose(sum(x[:-1]), x[-1], atol=0.001)
        if any(list(~sum_(np.array(list_)))):
            raise Exception('check_table_sums: ' + self.sample + ' has a mismatch in the total sums')

    def check_removals(self):
        """
        Check that the number of removed duplicates matches the number in the table and that the text of the second last
        line of the log are as expected
        """
        nums = list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", self.log_file[-2])))
        text = re.sub(r"\d+", "", self.log_file[-2])[:-1]
        # noinspection PyTypeChecker
        if ((text != 'samblaster: Removed         of      (.%) total read ids as duplicates using k memory in .S CPU seconds and S wall time.') |
            (int(nums[0]) != self.dups) |
            (any(int(i) < 0 for i in nums))):
            raise Exception('check_removals: ' + self.sample + ' has some issue (text or numeric related)')


class BaseRecalibrator(LogMain):
    """
    This class will check the baserecalibrator log
    """

    def __init__(self, path, sample):
        self.log_file = None
        self.log_template = None
        self.path = path
        self.sample = sample
        self.paired = self.single_paired()
        self.read_log()
        self._read_template()
        self.global_flags = []
        self.baserecalibrator = []  # Contains info on BaseRecalibrationEngine and BaseRecalibrator
        self.featuremanager = []
        self.progressmeter = []
        self.final_section = None
        self.split_log()

    def single_paired(self, table_path='data/fastq.csv'):
        """
        Check whether the sample is paired (R1 and R2) or single (only R1). This method is relevant
        for cases in which two different log files are generated
        :param table_path: Path in which we can find the fastq.csv (file containing this information)
        :return: Boolean value which will be stored as part of the class variables
        """

        df = pd.read_csv(table_path)

        bool_ = len(df[df.Sample == self.sample]) == 2

        return bool_

    def read_log(self):
        """
        Method to store the log file as part of the class variables
        """
        with open(self.path + self.sample + '_sort_nodup.recaldat.log') as f:
            self.log_file = f.readlines()

    def _read_template(self, path='data/template_recaldat.log'):
        """
        We read the global flag template so that we can compare this part more easily
        """
        with open(path) as f:
            self.log_template = f.readlines()

    def split_log(self):
        """
        This method is used to split the log into smaller and more manageable files (due to the large size of the file
        its better to split before to simplify our life). The outputs are stored as part of the class variables

        - ``self.global_flags``: Everything after [Global Flags] until that part ends
        - ``self.baserecalibrator``: Information on BaseRecalibrationEngine and BaseRecalibrator
        - ``self.featuremanager``: Information on FeatureManager
        - ``self.progressmeter``: Information on the ProgressMeter
        - ``self.final_section``: Last part of the log with final results
        """

        global_flags_bool = False
        baserecalibrator_bool = False
        featuremanager_bool = False
        progressmeter_bool = False

        for row in self.log_file:
            if row[:-1] == '[Global flags]':
                global_flags_bool = True
                continue

            if bool(re.search(r'(\d+-\d+-\d+)', row[:10])):
                global_flags_bool = False

            if bool(re.search(r'INFO  BaseRecalibrat', row)):
                baserecalibrator_bool = True

            if bool(re.search(r'read\(s\) filtered by:', row)):
                baserecalibrator_bool = True

            if bool(re.search(r'INFO  ProgressMeter', row)):
                progressmeter_bool = True

            if bool(re.search(r'INFO  FeatureManager', row)):
                featuremanager_bool = True

            if global_flags_bool:
                self.global_flags.append(row)
            elif baserecalibrator_bool:
                self.baserecalibrator.append(row)
                baserecalibrator_bool = False
            elif featuremanager_bool:
                self.featuremanager.append(row)
                featuremanager_bool = False
            elif progressmeter_bool:
                self.progressmeter.append(row)
                progressmeter_bool = False

        self.final_section = self.log_file[-11:]

    def check_log(self):
        """
        This method will run all the methods implemented for this class

        - ``check_running()``
        - ``check_correct_sample()``
        - ``check_global_flags_start()``
        - ``check_final_section()``
        - ``check_global_flags()``
        - ``check_baserecalibrator()``
        - ``check_featuremanager()``
        - ``check_progressmeter()``
        """

        self.check_running()
        self.check_correct_sample()
        self.check_global_flags_start()
        self.check_final_section()
        self.check_global_flags()
        self.check_baserecalibrator()
        self.check_featuremanager()
        self.check_progressmeter()

    def check_running(self):
        """
        Check that the log file contains a row written running:

        - ``Running:``
        """
        if self.log_file[1][:-1] != 'Running:':
            raise Exception('check_running: ' + self.sample + ' did not start running...')

    def check_correct_sample(self):
        """
        Make sure that the file that is being processed is the actual sample

        - If we are processing sample HSRR062650 we should only have this value in the command line string
        """
        if len(re.findall(self.sample, self.log_file[2][:-1])) == 0:
            raise Exception('check_correct_sample: ' + self.sample + ' should be processed however another sample has '
                                                                     'been processed instead')

    def check_global_flags_start(self):
        """
        Check we have a section with global flags:

        - ``[Global flags]``
        """
        if self.log_file[3][:-1] != '[Global flags]':
            raise Exception('check_global_flags: ' + self.sample + ' defining the global flags did not seem to work')

    def check_final_section(self):
        """
        Method to run checks on the final section part of the log, including the following methods:

        - ``check_final_section_length``
        - ``check_final_section_success``
        - ``check_final_section_others``

        """

        self.check_final_section_length()
        self.check_final_section_success()
        self.check_final_section_others()

    def check_global_flags(self):
        """
        Method to run checks on the global flags part of the log, including the following methods:

        - ``check_global_flags_length``
        - ``check_global_flags_variables``
        """

        self.check_global_flags_length()
        self.check_global_flags_variables()

    def check_baserecalibrator(self):
        """
        Method to run checks on the base recalibrator part of the log, including the following methods:

        - ``check_baserecalibrator_engine``
        - ``check_baserecalibrator_covariates``
        - ``check_baserecalibrator_filters``
        - ``check_baserecalibrator_quantization``
        """

        self.check_baserecalibrator_engine()
        self.check_baserecalibrator_covariates()
        self.check_baserecalibrator_filters()
        self.check_baserecalibrator_quantization()


    def check_featuremanager(self):
        """
        Method to run checks on the feature manager part of the log, including the following methods:

        - ``check_featuremanager_files``
        """

        self.check_featuremanager_files()



    def check_progressmeter(self):
        """
        Method to run checks on the progress meter part of the log, including the following methods:

        - ``check_progressmeter_len``
        - ``check_progressmeter_start_end``
        """

        self.check_progressmeter_len()
        self.check_progressmeter_start_end()

    def check_final_section_length(self):
        """
        Method to check the length of the final section

        - The length should be 11
        """
        if len(self.final_section) != 11:
            raise Exception('check_final_section_length: ' + self.sample + ' does not have the expected number of rows')

    def check_final_section_success(self):
        """
        Makes sure the Recalibration has been success
        """
        if self.final_section[1][:-1] != 'SUCCESS':
            raise Exception('check_final_section_success: ' + self.sample + ' has not been successful during '
                                                                            'BaseRecalibration')

    def check_final_section_others(self):
        """
        Check that the other rows also have the expected output, in particular we check that the following strings are
        present:

        - ``PSYoungGen``
        - ``ParOldGen``
        -  ``Metaspace``
        """
        s = ''.join(self.final_section[3:])
        if len(re.findall(r"\bPSYoungGen\b|\bParOldGen\b|\bMetaspace\b", s)) != 3:
            raise Exception('check_final_section_others: ' + self.sample + ' not all the expected fields are present '
                                                                           'in the last part of the log')

    def check_global_flags_length(self):
        """
        Method to check the length of the global flags section

        - The length should be 722
        """
        if len(self.global_flags) != 722:
            raise Exception('check_global_flags_length: ' + self.sample + ' does not have the expected number of rows')

    def check_global_flags_variables(self):
        """
        Method to check that all the global flags are set as expected

        - We compare the global flags with that of the template (check template_recaldat.log file for more details)
        - There are only 5 rows which are not supposed to be equal (they changed based on the computer session used)
            - ``uintx InitialHeapSize``
            - ``uintx MaxHeapSize``
            - ``uintx MaxNewSize``
            - ``uintx NewSize``
            - ``uintx OldSize``
        - All other variables should be the same
        """
        for row in range(len(self.global_flags)):
            original = self.global_flags[row]
            template = self.log_template[row]

            if (original != template) & (row not in [257, 304, 316, 349, 360]):
                raise Exception('check_global_flags_variables: ' + self.sample + ' does not have the right global '
                                                                                 'flags')

    def check_progressmeter_len(self):
        """
        Check correct length of this section:
        - For paired it should be 19 rows
        - For single it should be 18 rows
        """
        if self.paired:
            if len(self.progressmeter) != 19:
                raise Exception('check_progressmeter_len: ' + self.sample + ' does not have the right number of rows')
        else:
            if len(self.progressmeter) != 18:
                raise Exception('check_progressmeter_len: ' + self.sample + ' does not have the right number of rows')

    def check_progressmeter_start_end(self):
        """
        Check correct start and end statements
        - ``Starting traversal``
        - ``Traversal complete. Processed 180757 total reads in 2.6 minutes.``
        - We check that the strings are correct and also that the numeric part is larger than 0
        """
        start_text = self.progressmeter[0][-19:-1]
        end_text = re.findall('Traversal complete. Processed', self.progressmeter[-1][15:])
        end_nums = list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", self.progressmeter[-1][15:])))
        if ((start_text != 'Starting traversal') |
            (len(end_text) != 1) |
            (any(i <=0 for i in end_nums))):
            raise Exception('check_progressmeter_start_end: ' + self.sample + ' does not have the correct start and '
                                                                              'end statements in the ProgressMeter '
                                                                              'section')

    def check_featuremanager_files(self):
        """
        Check that the files used to extract the features are the correct ones. The files which should be used are:

        - ``hg38_resources/dbsnp_reannotated.vcf``
        - ``hg38_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf``
        - ``hg38_resources/1000G_omni2.5.hg38.vcf``
        - ``hg38_resources/wgs_calling_regions.hg38.interval_list``
        """

        list_ = ['hg38_resources/dbsnp_reannotated.vcf',
                 'hg38_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf',
                 'hg38_resources/1000G_omni2.5.hg38.vcf',
                 'hg38_resources/wgs_calling_regions.hg38.interval_list']

        s_check = ''.join(['\\' + 'b' + i + '\\' + 'b|' for i in list_])[:-1]
        s = ''.join(self.featuremanager)
        if len(re.findall(s_check, s)) != 4:
            raise Exception('check_featuremanager_files: ' + self.sample + ' did not get the features from the correct '
                                                                           'input files')

    def check_baserecalibrator_engine(self):
        """
        Method to check that the base recalibrator engine starts and ends correctly. We make sure that the following
        are present:

        - ``Initializing engine``
        - ``Done initializing engine``
        - ``Shutting down engine``
        """
        t1 = self.baserecalibrator[19][-20:-1] != 'Initializing engine'
        t2 = self.baserecalibrator[20][-25:-1] != 'Done initializing engine'
        t3 = self.baserecalibrator[-1][-21:-1] != 'Shutting down engine'

        if t1 | t2 | t3:
            raise Exception('check_baserecalibrator_engine: ' + self.sample + ' baserecalibrator engine did not work '
                                                                              'properly')

    def check_baserecalibrator_covariates(self):
        """
        Method to check that the base recalibrator engine considers the correct covariate values:

        - ``ReadGroupCovariate``
        - ``QualityScoreCovariate``
        - ``ContextCovariate``
        - ``CycleCovariate``
        """
        t1 = self.baserecalibrator[22][-19:-1] != 'ReadGroupCovariate'
        t2 = self.baserecalibrator[23][-22:-1] != 'QualityScoreCovariate'
        t3 = self.baserecalibrator[24][-17:-1] != 'ContextCovariate'
        t4 = self.baserecalibrator[25][-15:-1] != 'CycleCovariate'

        if t1 | t2 | t3 | t4:
            raise Exception('check_baserecalibrator_engine: ' + self.sample + ' not all the covariates have been used')

    def check_baserecalibrator_filters(self):
        """
        Method to check that the base recalibrator engine has execited the correct filters:

        - ``MappingQualityAvailableReadFilter``
        - ``MappedReadFilter``
        - ``NotSecondaryAlignmentReadFilter``
        - ``NotDuplicateReadFilter``
        - ``PassesVendorQualityCheckReadFilter``
        - ``WellformedReadFilter``
        """
        t1 = self.baserecalibrator[27][-35:-2] != 'MappingQualityAvailableReadFilter'
        t2 = self.baserecalibrator[28][-18:-2] != 'MappedReadFilter'
        t3 = self.baserecalibrator[29][-33:-2] != 'NotSecondaryAlignmentReadFilter'
        t4 = self.baserecalibrator[30][-24:-2] != 'NotDuplicateReadFilter'
        t5 = self.baserecalibrator[31][-36:-2] != 'PassesVendorQualityCheckReadFilter'
        t6 = self.baserecalibrator[32][-22:-2] != 'WellformedReadFilter'

        if t1 | t2 | t3 | t4 | t5 | t6:
            raise Exception('check_baserecalibrator_filters: ' + self.sample + ' not all the correct filters have been '
                                                                               'used on the base recalibrator')

    def check_baserecalibrator_quantization(self):
        """
        Method to check that the base recalibrator engine generates the report well and that the quantization passes

        - ``Calculating quantized quality scores...``
        - ``Writing recalibration report...``
        - ``...done!``
        - ``BaseRecalibrator was able to recalibrate 371510 reads``
        - We also check that the recalibration is done over a number greater than 0
        """
        t1 = self.baserecalibrator[33][-40:-1] != 'Calculating quantized quality scores...'
        t2 = self.baserecalibrator[34][-32:-1] != 'Writing recalibration report...'
        t3 = self.baserecalibrator[35][-9:-1] != '...done!'
        t4 = re.sub(r'[0-9]', '', self.baserecalibrator[36][-54:-1]) != 'BaseRecalibrator was able to recalibrate  reads'
        t5 = int(re.findall(r'\d+', self.baserecalibrator[36][-54:-1])[0]) < 0

        if t1 | t2 | t3 | t4 | t5:
            raise Exception('check_baserecalibrator_quantization: ' + self.sample + ' the quantization part did not '
                                                                                    'work properly')


class ApplyBQSR(LogMain):
    """
    This class will check the Apply BQSR log
    """

    def __init__(self, path, sample):
        self.log_file = None
        self.log_template = None
        self.path = path
        self.sample = sample
        self.paired = self.single_paired()
        self.read_log()
        self._read_template()
        self.global_flags = []
        self.applybqsr = []  # Contains info on BaseRecalibrationEngine and BaseRecalibrator
        self.featuremanager = []
        self.progressmeter = []
        self.final_section = None
        self.split_log()

    def single_paired(self, table_path='data/fastq.csv'):
        """
        Check whether the sample is paired (R1 and R2) or single (only R1). This method is relevant
        for cases in which two different log files are generated
        :param table_path: Path in which we can find the fastq.csv (file containing this information)
        :return: Boolean value which will be stored as part of the class variables
        """

        df = pd.read_csv(table_path)

        bool_ = len(df[df.Sample == self.sample]) == 2

        return bool_

    def read_log(self):
        """
        Method to store the log file as part of the class variables
        """
        with open(self.path + self.sample + '_sort_nodup.bqsr.log') as f:
            self.log_file = f.readlines()

    def _read_template(self, path='data/template_recaldat.log'):
        """
        We read the global flag template so that we can compare this part more easily
        """
        with open(path) as f:
            self.log_template = f.readlines()

    def split_log(self):
        """
        This method is used to split the log into smaller and more manageable files (due to the large size of the file
        its better to split before to simplify our life). The outputs are stored as part of the class variables

        - ``self.global_flags``: Everything after [Global Flags] until that part ends
        - ``self.applybqsr``: Information on ApplyBQSR and BaseRecalibrator
        - ``self.featuremanager``: Information on FeatureManager
        - ``self.progressmeter``: Information on the ProgressMeter
        - ``self.final_section``: Last part of the log with final results
        """

        global_flags_bool = False
        applybqsr_bool = False
        featuremanager_bool = False
        progressmeter_bool = False

        for row in self.log_file:
            if row[:-1] == '[Global flags]':
                global_flags_bool = True
                continue

            if bool(re.search(r'(\d+-\d+-\d+)', row[:10])):
                global_flags_bool = False

            if bool(re.search(r'INFO  ApplyBQSR', row)):
                applybqsr_bool = True

            if bool(re.search(r'INFO  ProgressMeter', row)):
                progressmeter_bool = True

            if bool(re.search(r'INFO  FeatureManager', row)):
                featuremanager_bool = True

            if global_flags_bool:
                self.global_flags.append(row)
            elif applybqsr_bool:
                self.applybqsr.append(row)
                applybqsr_bool = False
            elif featuremanager_bool:
                self.featuremanager.append(row)
                featuremanager_bool = False
            elif progressmeter_bool:
                self.progressmeter.append(row)
                progressmeter_bool = False

        self.final_section = self.log_file[-9:]

    def check_log(self):
        """
        This method will run all the methods implemented for this class

        - ``check_running()``
        - ``check_correct_sample()``
        - ``check_global_flags_start()``
        - ``check_final_section()``
        - ``check_global_flags()``
        - ``check_applybqsr()``
        - ``check_featuremanager()``
        - ``check_progressmeter()``
        """

        self.check_running()
        self.check_correct_sample()
        self.check_global_flags_start()
        self.check_final_section()
        self.check_global_flags()
        self.check_applybqsr()
        self.check_featuremanager()
        self.check_progressmeter()

    def check_running(self):
        """
        Check that the log file contains a row written running:

        - ``Running:``
        """
        if self.log_file[1][:-1] != 'Running:':
            raise Exception('check_running: ' + self.sample + ' did not start running...')

    def check_correct_sample(self):
        """
        Make sure that the file that is being processed is the actual sample

        - If we are processing sample HSRR062650 we should only have this value in the command line string
        """
        if len(re.findall(self.sample, self.log_file[2][:-1])) == 0:
            raise Exception('check_correct_sample: ' + self.sample + ' should be processed however another sample has '
                                                                     'been processed instead')

    def check_global_flags_start(self):
        """
        Check we have a section with global flags:

        - ``[Global flags]``
        """
        if self.log_file[3][:-1] != '[Global flags]':
            raise Exception('check_global_flags: ' + self.sample + ' defining the global flags did not seem to work')

    def check_final_section(self):
        """
        Method to run checks on the final section part of the log, including the following methods:

        - ``check_final_section_length``
        - ``check_final_section_success``
        - ``check_final_section_others``

        """
        # DONE
        self.check_final_section_length()
        self.check_final_section_others()

    def check_global_flags(self):
        """
        Method to run checks on the global flags part of the log, including the following methods:

        - ``check_global_flags_length``
        - ``check_global_flags_variables``
        """
        # DONE
        self.check_global_flags_length()
        self.check_global_flags_variables()

    def check_applybqsr(self):
        """
        Method to run checks on the ApplyBQSR part of the log, including the following methods:

        - ``check_applybqsr_engine``
        - ``check_applybqsr_quantization``
        """

        self.check_applybqsr_engine()
        self.check_applybqsr_quantization()

    def check_featuremanager(self):
        """
        Method to run checks on the feature manager part of the log, including the following methods:

        - ``check_featuremanager_files``
        """

        self.check_featuremanager_files()

    def check_progressmeter(self):
        """
        Method to run checks on the progress meter part of the log, including the following methods:

        - ``check_progressmeter_len``
        - ``check_progressmeter_start_end``
        """

        self.check_progressmeter_len()
        self.check_progressmeter_start_end()

    def check_final_section_length(self):
        """
        Method to check the length of the final section

        - The length should be 11
        """
        if len(self.final_section) != 9:
            raise Exception('check_final_section_length: ' + self.sample + ' does not have the expected number of rows')

    def check_final_section_others(self):
        """
        Check that the other rows also have the expected output, in particular we check that the following strings are
        present:

        - ``PSYoungGen``
        - ``ParOldGen``
        - ``Metaspace``
        """
        s = ''.join(self.final_section[1:])
        if len(re.findall(r"\bPSYoungGen\b|\bParOldGen\b|\bMetaspace\b", s)) != 3:
            raise Exception('check_final_section_others: ' + self.sample + ' not all the expected fields are present '
                                                                           'in the last part of the log')

    def check_global_flags_length(self):
        """
        Method to check the length of the global flags section

        - The length should be 722
        """
        if len(self.global_flags) != 722:
            raise Exception('check_global_flags_length: ' + self.sample + ' does not have the expected number of rows')

    def check_global_flags_variables(self):
        """
        Method to check that all the global flags are set as expected

        - We compare the global flags with that of the template (check template_recaldat.log file for more details)
        - There are only 5 rows which are not supposed to be equal (they changed based on the computer session used)
            - ``uintx InitialHeapSize``
            - ``uintx MaxHeapSize``
            - ``uintx MaxNewSize``
            - ``uintx NewSize``
            - ``uintx OldSize``
        - All other variables should be the same
        """
        for row in range(len(self.global_flags)):
            original = self.global_flags[row]
            template = self.log_template[row]

            if (original != template) & (row not in [257, 304, 316, 349, 360]):
                raise Exception('check_global_flags_variables: ' + self.sample + ' does not have the right global '
                                                                                 'flags')

    def check_progressmeter_len(self):
        """
        Check correct length of this section:

        - For paired and single it should be 4 rows

        """
        if len(self.progressmeter) != 4:
            raise Exception('check_progressmeter_len: ' + self.sample + ' does not have the right number of rows')

    def check_progressmeter_start_end(self):
        """
        Check correct start and end statements

        - ``Starting traversal``
        - ``Traversal complete. Processed 180757 total reads in 2.6 minutes.``
        - We check that the strings are correct and also that the numeric part is larger than 0
        """
        start_text = self.progressmeter[0][-19:-1]
        end_text = re.findall('Traversal complete. Processed', self.progressmeter[-1][15:])
        end_nums = list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", self.progressmeter[-1][15:])))
        if ((start_text != 'Starting traversal') |
            (len(end_text) != 1) |
            (any(i <=0 for i in end_nums))):
            raise Exception('check_progressmeter_start_end: ' + self.sample + ' does not have the correct start and '
                                                                              'end statements in the ProgressMeter '
                                                                              'section')

    def check_featuremanager_files(self):
        """
        Check that the files used to extract the features are the correct ones. The files which should be used are:

        - ``hg38_resources/wgs_calling_regions.hg38.interval_list``
        """

        list_ = ['hg38_resources/wgs_calling_regions.hg38.interval_list']

        s_check = ''.join(['\\' + 'b' + i + '\\' + 'b|' for i in list_])[:-1]
        s = ''.join(self.featuremanager)
        if len(re.findall(s_check, s)) != len(list_):
            raise Exception('check_featuremanager_files: ' + self.sample + ' did not get the features from the correct '
                                                                           'input files')

    def check_applybqsr_engine(self):
        """
        Method to check that the applybqsr engine starts and ends correctly. We make sure that the following
        are present:

        - ``Initializing engine``
        - ``Done initializing engine``
        - ``Shutting down engine``
        """
        t1 = self.applybqsr[19][-20:-1] != 'Initializing engine'
        t2 = self.applybqsr[20][-25:-1] != 'Done initializing engine'
        t3 = self.applybqsr[-1][-21:-1] != 'Shutting down engine'

        if t1 | t2 | t3:
            raise Exception('check_applybqsr_engine: ' + self.sample + ' applybqsr engine did not work '
                                                                              'properly')

    def check_applybqsr_quantization(self):
        """
        Method to check that the applybqsr engine generates the report well and that the quantization passes
        """
        t1 = self.applybqsr[-2][-22:-2] != 'WellformedReadFilter'

        if t1:
            raise Exception('check_applybqsr_quantization: ' + self.sample + ' the quantization part did not '
                                                                                    'work properly')








