import os
import argparse
import configparser

class ArgumentParser:
    def __init__(self):
        self.parse_program_arguments()
        self.config = configparser.ConfigParser()
        self.base_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        self.config_ini = os.path.join(self.base_dir, 'settings', 'config.ini')
        self.config.read(self.config_ini)

        if self.args.command == 'settings':
            print('Applying settings...')
            self.apply_settings_to_config(args, config)

        if self.args.automatic:
            self.apply_defaults_from_config('VCF')
            self.apply_defaults_from_config('BSA')
        elif 'vcf_generator' in self.args:
            self.apply_defaults_from_config('VCF')
        elif 'analysis' in self.args:
            self.apply_defaults_from_config('BSA')
 
    def apply_settings_to_config(self, args, config):
        for arg in vars(args):
            if arg.startswith('set_') and getattr(args, arg) is not None:
                value = getattr(args, arg)
                config_key = arg.replace('set_', '')
                for section in config.sections():
                    if config.has_option(section, config_key):
                        print(f'Default {config_key} in section {section} is now set to {value}')
                        config.set(section, config_key, str(value))
                        break  # break after finding the option in the section
        with open('config.ini', 'w') as configfile:
            config.write(configfile)
        quit()

    def apply_defaults_from_config(self, section):
        for arg in vars(self.args).keys():
            if getattr(self.args, arg) is None and self.config.has_option(section, arg):
                value = self.config.get(section, arg)
                print(f'Default applied: {arg}:{value}')
                setattr(self.args, arg, value)

    def add_bsa_arguments(self, parser):
        bsa_options = parser.add_argument_group('BSA analysis options', 'Options for BSA analysis. Defaults can be changed using the settings positional argument. phytobsa settings -h for for info')
        bsa_options.add_argument('-ls', '--loess_span', type=float, default=None, help="Influences smoothing parameters.")
        bsa_options.add_argument('-si', '--shuffle_iterations', type=int, default=None, help="Iterations of bootstrapping during empirical cutoff calculations. Below 1000 can yield inconsistent results")
        bsa_options.add_argument('-sb', '--smooth_edges_bounds', type=int, default=None, help="Number of mirrored datapoints at chromosome edges to correct for loess edge bias. Increase if edge bias seems high")
        bsa_options.add_argument('-fin', '--filter_indels', default=None, type=str, help="Filter out insertion-deletion mutations.")
        bsa_options.add_argument('-fems', '--filter_ems', default=None, type=str, help="Filter results to only include mutations likely to arise from EMS treatment")
        bsa_options.add_argument('-snpmsk', '--snpmask_path', default=None, type=str, help="Path to VCF file containing background snps.")

    def add_vcf_gen_arguments(self, parser):
        vcf_gen_options = parser.add_argument_group('VCF generation options', 'Options for VCF generation. Defaults can be changed using the settings positional argument phytobsa settings -h for more info')
        vcf_gen_options.add_argument('-rgn', '--reference_genome_name', required=False, default=None, type=str, help='Reference genome name')
        vcf_gen_options.add_argument('-ssdb', '--snpEff_species_db', default=None, type=str, help = 'The name of your snpEff database name.')
        vcf_gen_options.add_argument('-rgs', '--reference_genome_source', default=None, type=str, help = 'Optional, if you wish the pipeline to download your reference from a url')
        vcf_gen_options.add_argument('-ks','--known_snps', default=None, type=str, help = 'VCF file containing background SNPs. Helps improve output quality')
        vcf_gen_options.add_argument('-t', '--threads_limit', default=None, type=str, help='Maximum threads you wish to use for analysis')
        vcf_gen_options.add_argument('-p', '--call_variants_in_parallel', default=None, type=bool, help='Run gatk haplotype caller in parallel')
        vcf_gen_options.add_argument('-c','--cleanup', default=None, type=bool, help='If true, intermediate files will de deleted. False for troubleshooting and archiving files.' )

    def parse_program_arguments(self):
        
        # Main parser
        main_parser = argparse.ArgumentParser(description='PyAtBSA main command line script...')
        main_parser.add_argument('-a', '--automatic', action='store_true', help='Run PhytoBSA in automatic mode. Primary, default workflow.')
        self.add_bsa_arguments(main_parser)
        self.add_vcf_gen_arguments(main_parser)

        subparsers = main_parser.add_subparsers(help='Sub-command help', dest='command')
        ## Analysis subparser
        parser_analysis = subparsers.add_parser('analysis', help='Run PhytoBSA analysis seperately.')
        ### REQUIRED
        parser_analysis.add_argument('-n', '--name', required=True, type=str, help='name of the line you wish to analyze. Will be used to name output files.')
        parser_analysis.add_argument('-vt', '--vcf_table_path', required=True, type=str, help='path to the vcf table you wish to analyze.')
        parser_analysis.add_argument('-st', '--segregation_type', required=True, type=str, help="Recessive(R) or Dominant(D)")
        self.add_bsa_arguments(parser_analysis)

        ## VCF generation subparser
        parser_vcf_gen = subparsers.add_parser('vcf_generator', help='Run PhytoBSA VCF Generator seperately')
        ### REQUIRED
        parser_vcf_gen.add_argument('-wt', '--wt_input', required=True, default=None, type=str, help='Wild-type bulk fasta file(s)')
        parser_vcf_gen.add_argument('-mu', '--mu_input', required=True, default=None, type=str, help='mutant bulk fast files')
        ### OPTIONAL - the below arguments can be input via command line or sourced from ./settings/config.ini
        self.add_vcf_gen_arguments(parser_vcf_gen)

        ## Log database subparser
        parser_log_db = subparsers.add_parser('logdb', help='Gather log information from log database')
        parser_log_db.add_argument('-vcf', "--vcf_ulid_log", default = None, type = str, help="Retrieve vcf log information based on vcf ulid input")
        parser_log_db.add_argument('-an', "--analysis_ulid_log", default = None, type = str, help="Retrieve analysis log information based on analysis ulid input")
        parser_log_db.add_argument('-name', "--line_name_log", default = None, type = str, help="Retrieve all information associated with the provided line name. ")    
        parser_log_db.add_argument('-core', '--core_ulid_log', default = None, type = str, help="Retrieve all information associated with the provided core log ulid")
        
        ## Settings subparser
        parser_settings = subparsers.add_parser('settings', help='Update default settings.')
        parser_settings.add_argument('--set_data_dir', default=None, type=str, help='set Data directory. This must be set for program to run')
        #vcf_gen default run settings. 
        vcf_settings = parser_settings.add_argument_group('VCF generation default settings' 'These settings will be automatically applied if not explicity provided in automatic or VCF generation mode')
        vcf_settings.add_argument('--set_reference_genome_name', required=False, default=None, type=str, help='Set default reference genome name')
        vcf_settings.add_argument('--set_snpEff_species_db', default=None, type=str, help = 'Set default snpEff database name.')
        vcf_settings.add_argument('--set_reference_genome_source', default=None, type=str, help = 'Set default reference genome source')
        vcf_settings.add_argument('--set_known_snps', default=None, type=str, help = 'Set default VCF file containing background SNPs.')
        vcf_settings.add_argument('--set_threads_limit', default=None, type=str, help='Set default maximum threads for analysis')
        vcf_settings.add_argument('--set_call_variants_in_parallel', default=None, type=bool, help='Set default for running gatk haplotype caller in parallel')
        vcf_settings.add_argument('--set_cleanup', default=None, type=bool, help='Set default for cleanup. If true, intermediate files will de deleted. False for troubleshooting and archiving files.' )
        #BSA default run settings.
        bsa_settings = parser_settings.add_argument_group('BSA default settings', 'These settings will be automatically applied if not explicitly passed to automatic or BSA mode.')
        bsa_settings.add_argument('--set_loess_span', type=float, help="Set default loess_span.")
        bsa_settings.add_argument('--set_shuffle_iterations', type=int, help="Set default shuffle_iterations.")
        bsa_settings.add_argument('--set_smooth_edges_bounds', type=int, help="Set default smooth_edges_bounds.")
        bsa_settings.add_argument('--set_filter_indels', type=bool, help="Set default filter_indels.")
        bsa_settings.add_argument('--set_filter_ems', type=bool, help="Set default filter_ems.")
        bsa_settings.add_argument('--set_snpmask_path', type=str, help="Set default snp mask VCF path. VCF should contain known background SNPs")
        
        #parse args
        self.args = main_parser.parse_args()
