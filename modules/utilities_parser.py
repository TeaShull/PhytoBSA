import os
import argparse
import configparser
import ast


class ArgumentParser:
    
    def __init__(self):
        self.parse_program_arguments()
        self.config = configparser.ConfigParser()
        self.base_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        self.config_ini = os.path.join(self.base_dir, 'settings', 'config.ini')
        self.config.read(self.config_ini)

        if self.args.command == 'settings':
            if not any(vars(self.args).values()):
                parser.print_help()
            elif self.args.list:
                self.print_settings_from_config()
            else:
                print('Applying settings...')
                self.apply_settings_to_config()
        
        if self.args.automatic:
            self.apply_defaults_from_config('GENERAL')
            self.apply_defaults_from_config('VCF')
            self.apply_defaults_from_config('BSA')
        
        elif self.args.command == 'vcf_generator':
            self.apply_defaults_from_config('GENERAL')
            self.apply_defaults_from_config('VCF')
        
        elif self.args.command == 'analysis':
            self.apply_defaults_from_config('GENERAL')
            self.apply_defaults_from_config('BSA')

    def print_settings_from_config(self):
        for section in self.config.sections():
            print(f'[{section}]')
            for key, value in self.config.items(section):
                print(f'{key} = {value}')
            print(" ")
        quit()

    def apply_settings_to_config(self):
        for arg in vars(self.args):
            if arg.startswith('set_') and getattr(self.args, arg) is not None:
                value = str(getattr(self.args, arg))
                config_key = arg.replace('set_', '')
                print(f'{value}:{config_key}')
                for section in self.config.sections():
                    if self.config.has_option(section, config_key):
                        print(f'Default updated |{section}|{config_key} = {value}')
                        self.config.set(section, config_key, value)
                        break
        with open(self.config_ini, 'w') as configfile:
            self.config.write(configfile)
        quit()

    def apply_defaults_from_config(self, section):
        for arg in vars(self.args).keys():
            if getattr(self.args, arg) is None and self.config.has_option(section, arg):
                value = self.config.get(section, arg)
                try:
                    value = ast.literal_eval(value)
                except (ValueError, SyntaxError):
                    pass
                print(f'Default applied: {arg}:{value}')
                setattr(self.args, arg, value)

    def add_bsa_arguments(self, parser):
        bsa_options = parser.add_argument_group('BSA analysis options', 'Options for BSA analysis. Defaults can be changed using the settings positional argument. phytobsa settings -h for for info')
        bsa_options.add_argument('-ls', '--loess_span', type=float, default=None, help="Influences smoothing parameters.")
        bsa_options.add_argument('-si', '--shuffle_iterations', type=int, default=None, help="Iterations of bootstrapping during empirical cutoff calculations. Below 1000 can yield inconsistent results")
        bsa_options.add_argument('-sb', '--smooth_edges_bounds', type=int, default=None, help="Number of mirrored datapoints at chromosome edges to correct for loess edge bias. Increase if edge bias seems high")
        bsa_options.add_argument('-fin', '--filter_indels', default=None, type=str, help="Filter out insertion-deletion mutations.")
        bsa_options.add_argument('-fems', '--filter_ems', default=None, type=str, help="Filter results to only include mutations likely to arise from EMS treatment")
        bsa_options.add_argument('-rco', '--ratio_cutoff', default=None, type=float, help="Used to filter results based on a ratio cutoff number. Increase to 0.2 or 0.3 if there is alot of noise at lower ratio bounds")
        bsa_options.add_argument('-msk', '--mask_snps', default=None, type=bool, help="Set to true if you have a snpmask file configured and would like to mask known snps in your analysis.")
        bsa_options.add_argument('-cc', '--critical_cutoff', default=None, type=float, help="Set the critical cutoff value for what is considered a significant polymorphism.")
    
    def add_vcf_gen_arguments(self, parser):
        vcf_gen_options = parser.add_argument_group('VCF generation options', 'Options for VCF generation. Defaults can be changed using the settings positional argument phytobsa settings -h for more info')
        vcf_gen_options.add_argument('-p', '--call_variants_in_parallel', default=None, type=bool, help='Run gatk haplotype caller in parallel')
        vcf_gen_options.add_argument('-c','--cleanup', default=None, type=bool, help='If true, intermediate files will de deleted. False for troubleshooting and archiving files.' )
        vcf_gen_options.add_argument('-cft', '--cleanup_filetypes', default=None, type=list, help="Filetypes to clean out after VCF generation is complete. format - ['*file_suffix', exc] example - ['*.tmp', '*.metrics']")
        vcf_gen_options.add_argument('-ocp', '--omit_chrs_patterns', default=None, type=list, help='Header patterns to omit from reference chromosomes. Useful for removing >mt(mitochondrial) and other unneeded reference sequences')

    def add_reference_name_argument(self, parser):
        parser.add_argument('-r', '--reference_name', required=False, type=str, help='Name of the reference genome')

    def parse_program_arguments(self):
        
        # Main parser
        main_parser = argparse.ArgumentParser(description='PyAtBSA main command line script...')
        main_parser.add_argument('-a', '--automatic', action='store_true', help='Run PhytoBSA in automatic mode. Primary, default workflow.')
        self.add_bsa_arguments(main_parser)
        self.add_vcf_gen_arguments(main_parser)
        self.add_reference_name_argument(main_parser)

        subparsers = main_parser.add_subparsers(help='Sub-command help', dest='command')
        ## Analysis subparser
        parser_analysis = subparsers.add_parser('analysis', help='Run PhytoBSA analysis seperately.', add_help=True)
        ### REQUIRED
        parser_analysis.add_argument('-n', '--name', required=True, type=str, help='name of the line you wish to analyze. Will be used to name output files.')
        parser_analysis.add_argument('-vt', '--vcf_table_path', required=True, type=str, help='path to the vcf table you wish to analyze.')
        ##Optional but recommended 
        parser_analysis.add_argument('-st', '--segregation_type', required=False, default=None, type=str, help="Recessive(R) or Dominant(D)")
        self.add_bsa_arguments(parser_analysis)
        self.add_reference_name_argument(parser_analysis)

        ## VCF generation subparser
        parser_vcf_gen = subparsers.add_parser('vcf_generator', help='Run PhytoBSA VCF Generator seperately',add_help=True)
        ### REQUIRED
        parser_vcf_gen.add_argument('-n', '--name', required=True, type=str, help='name of the line you are generating a VCF file for. Will be used to name output files.')
        parser_vcf_gen.add_argument('-wt', '--wt_input', required=True, default=None, type=str, help='Wild-type bulk fasta file(s)')
        parser_vcf_gen.add_argument('-mu', '--mu_input', required=True, default=None, type=str, help='mutant bulk fast files')
        ### OPTIONAL - the below arguments can be input via command line or sourced from ./settings/config.ini
        self.add_vcf_gen_arguments(parser_vcf_gen)
        self.add_reference_name_argument(parser_vcf_gen)

        ## Log database subparser
        parser_log_db = subparsers.add_parser('logdb', help='Gather log information from log database', add_help=True)
        parser_log_db.add_argument('-vcf', "--vcf_ulid_log", default = None, type = str, help="Retrieve vcf log information based on vcf ulid input")
        parser_log_db.add_argument('-an', "--analysis_ulid_log", default = None, type = str, help="Retrieve analysis log information based on analysis ulid input")
        parser_log_db.add_argument('-name', "--line_name_log", default = None, type = str, help="Retrieve all information associated with the provided line name. ")    
        parser_log_db.add_argument('-core', '--core_ulid_log', default = None, type = str, help="Retrieve all information associated with the provided core log ulid")
        
        ## Settings subparser
        parser_settings = subparsers.add_parser('settings', help='Update default settings.', add_help=True)
        parser_settings.add_argument('--set_data_dir', default=None, type=str, help='set Data directory. This must be set for program to run')
        parser_settings.add_argument('--set_threads_limit', default=None, type=int, help="Set the threads limit for BSA and for VCF generation. If not set, threads will be detected and threads -2 will be used. ")
        parser_settings.add_argument('--set_reference_name', default=None, required=False, type=str, help='Name of the reference genome')
        parser_settings.add_argument('--list', default=False, action='store_true', help='List settings')
        
        #vcf_gen default run settings. 
        vcf_settings = parser_settings.add_argument_group('VCF generation default settings' 'These settings will be automatically applied if not explicity provided in automatic or VCF generation mode')
        vcf_settings.add_argument('--set_call_variants_in_parallel', default=None, type=bool, help='Set default for running gatk haplotype caller in parallel')
        vcf_settings.add_argument('--set_cleanup', default=None, type=bool, help='Set default for cleanup. If true, intermediate files will de deleted. False for troubleshooting and archiving files.' )
        vcf_settings.add_argument('--set_cleanup_filetypes', default=None, type=list, help="set default for cleanup filetypes. ordered list of globs for files you wish to clear out after vcf generation process")
        vcf_settings.add_argument('--set_omit_chrs_patterns', default=None, type=list, help="set defaults for filtering reference chromosome contigs. Useful for filtering non-genomic reference contigs to speed up vcf generation")
        #BSA default run settings.
        bsa_settings = parser_settings.add_argument_group('BSA default settings', 'These settings will be automatically applied if not explicitly passed to automatic or BSA mode.')
        bsa_settings.add_argument('--set_loess_span', type=float, help="Set default loess_span.")
        bsa_settings.add_argument('--set_shuffle_iterations', type=int, help="Set default shuffle_iterations.")
        bsa_settings.add_argument('--set_smooth_edges_bounds', type=int, help="Set default smooth_edges_bounds.")
        bsa_settings.add_argument('--set_filter_indels', type=bool, help="Set default filter_indels.")
        bsa_settings.add_argument('--set_filter_ems', type=bool, help="Set default filter_ems.")
        bsa_settings.add_argument('--set_ratio_cutoff', type=float, help="Set default ratio cutoff bound.")
        bsa_settings.add_argument('--set_mask_snps', type=bool, help="set default mask_snps boolean value.")
        bsa_settings.add_argument('--set_critical_cutoff', type=float, help="Set default critical cutoff value.")

        #parse args
        self.args = main_parser.parse_args()
