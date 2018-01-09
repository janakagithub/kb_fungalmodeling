package kb_fungalmodeling::kb_fungalmodelingImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org
our $VERSION = '0.0.1';
our $GIT_URL = 'https://github.com/janakagithub/kb_fungalmodeling.git';
our $GIT_COMMIT_HASH = '661a64a87ffc7ed21054f1428623f1182b957f70';

=head1 NAME

kb_fungalmodeling

=head1 DESCRIPTION

A KBase module: kb_fungalmodeling
This module  build fungal models based on fungal genomes.

=cut

#BEGIN_HEADER
use Bio::KBase::AuthToken;
use Workspace::WorkspaceClient;
use AssemblyUtil::AssemblyUtilClient;
use KBaseReport::KBaseReportClient;
use GenomeProteomeComparison::GenomeProteomeComparisonClient;
use fba_tools::fba_toolsClient;
use Config::IniFiles;
use Bio::SeqIO;
use UUID::Random;
use Data::Dumper;
#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR

    my $config_file = $ENV{ KB_DEPLOYMENT_CONFIG };
    my $cfg = Config::IniFiles->new(-file=>$config_file);
    my $scratch = $cfg->val('kb_fungalmodeling', 'scratch');
    my $wsInstance = $cfg->val('kb_fungalmodeling','workspace-url');
    my $callbackURL = $ENV{ SDK_CALLBACK_URL };

    $self->{scratch} = $scratch;
    $self->{callbackURL} = $callbackURL;
    $self->{'workspace-url'} = $wsInstance;
    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 build_fungal_model

  $output = $obj->build_fungal_model($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_fungalmodeling.fungalmodelbuiltInput
$output is a kb_fungalmodeling.fungalmodelbuiltOutput
fungalmodelbuiltInput is a reference to a hash where the following keys are defined:
	workspace has a value which is a string
	genome_ref has a value which is a string
	template_model has a value which is a string
	gapfill_model has a value which is an int
	translation_policy has a value which is a string
	output_model has a value which is a string
fungalmodelbuiltOutput is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string

</pre>

=end html

=begin text

$params is a kb_fungalmodeling.fungalmodelbuiltInput
$output is a kb_fungalmodeling.fungalmodelbuiltOutput
fungalmodelbuiltInput is a reference to a hash where the following keys are defined:
	workspace has a value which is a string
	genome_ref has a value which is a string
	template_model has a value which is a string
	gapfill_model has a value which is an int
	translation_policy has a value which is a string
	output_model has a value which is a string
fungalmodelbuiltOutput is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string


=end text



=item Description



=back

=cut

sub build_fungal_model
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to build_fungal_model:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'build_fungal_model');
    }

    my $ctx = $kb_fungalmodeling::kb_fungalmodelingServer::CallContext;
    my($output);
    #BEGIN build_fungal_model

    print("Starting fungal model building method. Parameters:\n");
    print(Dumper($params) . "\n");


    my $fbaO = new fba_tools::fba_toolsClient( $self->{'callbackURL'},
                                                            ( 'service_version' => 'release',
                                                              'async_version' => 'release',
                                                            )
                                                           );

    my $protC = new GenomeProteomeComparison::GenomeProteomeComparisonClient( $self->{'callbackURL'},
                                                            ( 'service_version' => 'release',
                                                              'async_version' => 'release',
                                                            )
                                                           );

    my $reportHandle = new KBaseReport::KBaseReportClient( $self->{'callbackURL'},
                                                            ( 'service_version' => 'release',
                                                              'async_version' => 'release',
                                                            )
                                                          );

    my $token=$ctx->token;
    my $wshandle= Workspace::WorkspaceClient->new($self->{'workspace-url'},token=>$token);

    my $template_ws = 'jplfaria:narrative_1510597445008';# 'janakakbase:narrative_1513399583946'; # template workspaces
    my $template_genome_ref = 'fungal_template.genome';
    my $template_model_ref = 'master_fungal_template';
    my $ws_name = $params->{workspace};
    my $protCompId = 'proteinComp'.$params->{genome_ref};
    my $tmpGenome;
    my $tmpModel;

    my $templateId = {
      default_temp => [$template_model_ref, $template_genome_ref],
      iJL1454 => ['iJL1454_KBase', 'Aspergillus_terreus_NIH2624'],
      iNX804  => ['iNX804_KBase','Candida_glabrata_ASM254'],
      iCT646 => ['iCT646_KBase','Candida_tropicali_MYA-3404'],
      iOD907 => ['iOD907_KBase','Kluyveromyces_lactis_NRRL'],
      iJDZ836 => ['iJDZ836_KBase','Neurospora_crassa_OR74A'],
      Yeast => ['iMM904_KBase','Saccharomyces_cerevisiae_5288c']

    };

  my $eachTemplateHash;
  my $eachModelRxns;
  my $eachModelMSRxns;
  my $eachModelCpds;
  my $eachModelMSCpds;
  my @newModelArr;

# Generate stats
  foreach my $k (keys $templateId){
    if ($k eq 'default_temp'){
      next;
    }
    else{
        eval {
           print "retrieving individual models from template $k\n";
           my $eachTemplate = $wshandle->get_objects([{workspace=>$template_ws,name=>$templateId->{$k}->[0]}] )->[0]{data};#{modelreactions};

           for (my $i=0; $i< @{$eachTemplate->{modelreactions}}; $i++){

            $eachTemplateHash->{$eachTemplate->{modelreactions}->[$i]->{id}} = [$k, $templateId->{$k}->[1]];
            my @msr1 = split /\//, $eachTemplate->{modelreactions}->[$i]->{reaction_ref};
            my @msr = split /_/, $msr1[-1];

            if ($msr[0] eq 'rxn00000'){

              $eachModelRxns->{$templateId->{$k}->[1]}++;
            }
            else {

              $eachModelRxns->{$templateId->{$k}->[1]}++;
              $eachModelMSRxns->{$templateId->{$k}->[1]}++;
            }

           }

           for (my $i=0; $i< @{$eachTemplate->{modelcompounds}}; $i++){
            my @msc1 = split /\//, $eachTemplate->{modelcompounds}->[$i]->{compound_ref};
            my @msc = split /_/, $msc1[-1];

            if ($msc[0] eq 'cpd00000'){

              $eachModelCpds->{$templateId->{$k}->[1]}++;
            }
            else {

              $eachModelCpds->{$templateId->{$k}->[1]}++;
              $eachModelMSCpds->{$templateId->{$k}->[1]}++;
            }

           }

        };
        if ($@) {
           die "Error loading object from the workspace:\n".$@;
        }
    }

  }

  print &Dumper ($eachModelRxns);
  print &Dumper ($eachModelMSRxns);
  print &Dumper ($eachModelCpds);
  print &Dumper ($eachModelMSCpds);

    if (defined $params->{template_model}) {

      print "Template selected as $params->{template_model} \n";
      $tmpModel = $templateId->{$params->{template_model}}->[0];
      $tmpGenome = $templateId->{$params->{template_model}}->[1];
      print &Dumper ($templateId);

    }
    my $tran_policy;
    if (defined $params->{translation_policy}){
      $tran_policy = $params->{translation_policy};

    }

    print "producing a proteome comparison object between $params->{genome_ref} and $tmpGenome\n";

    my $protComp =  $protC->blast_proteomes({
        genome1ws => $params->{workspace},
        genome1id => $params->{genome_ref},
        genome2ws => $template_ws,
        genome2id => $tmpGenome,
        output_ws => $params->{workspace},
        output_id => $protCompId
    });

    print "Producing a draft model based on $protCompId proteome comparison\n";


    my $gpModelFromSource;
    if ($params->{gapfill_model} == 1){

        my $dr_model =$params->{genome_ref}."_draftModel";
        my $fba_modelProp =  $fbaO->propagate_model_to_new_genome({
            fbamodel_id => $tmpModel,
            fbamodel_workspace => $template_ws,
            proteincomparison_id => $protCompId,
            proteincomparison_workspace => $params->{workspace},
            fbamodel_output_id =>  $dr_model,
            workspace => $params->{workspace},
            keep_nogene_rxn => 1,
            #media_id =>
            #media_workspace =>
            minimum_target_flux => 0.1,
            translation_policy => $tran_policy
            #output_id =>  $dr_model
        });

        print "Running gapfill from the source model, may take a while....\n ";

        $gpModelFromSource = $fbaO->gapfill_metabolic_model ({

            fbamodel_id => $dr_model,
            fbamodel_output_id => $params->{output_model},
            source_fbamodel_id => $tmpModel,
            source_fbamodel_workspace => $template_ws,
            workspace => $params->{workspace},
            target_reaction => 'bio1',
            feature_ko_list => [],
            reaction_ko_list => [],
            custom_bound_list => [],
            media_supplement_list =>  [],
            minimum_target_flux => 0.1
        });
    }
    else {


        my $fba_modelProp =  $fbaO->propagate_model_to_new_genome({
            fbamodel_id => $tmpModel,
            fbamodel_workspace => $template_ws,
            proteincomparison_id => 'proteinCompAspergillus_oryzae_BCC7051', #$protCompId,
            proteincomparison_workspace => $params->{workspace},
            fbamodel_output_id =>  $params->{output_model},
            workspace => $params->{workspace},
            keep_nogene_rxn => 1,
            #media_id =>
            #media_workspace =>
            minimum_target_flux => 0.1,
            translation_policy => $tran_policy
            #output_id =>  $params->{output_model}
        });
    }

       print &Dumper ($gpModelFromSource);


    eval {
        #my $newModel = $wshandle->get_objects([{workspace=>'janakakbase:narrative_1509987427391',name=>'Psean1_DF_GP'}])->[0]{data}->{modelreactions};
        my $newModel = $wshandle->get_objects([{workspace=>$params->{workspace},name=>$params->{output_model}}])->[0]{data};# ->{modelreactions};
        for (my $i=0; $i< @{$newModel->{modelreactions}}; $i++){

            push (@newModelArr,$newModel->{modelreactions}->[$i]->{id} );

            my @msr1 = split /\//, $newModel->{modelreactions}->[$i]->{reaction_ref};
            my @msr = split /_/, $msr1[-1];

            if ($msr[0] eq 'rxn00000'){

              $eachModelRxns->{ $params->{output_model}}++;
            }
            else {

              $eachModelRxns->{ $params->{output_model}}++;
              $eachModelMSRxns->{ $params->{output_model}}++;
            }

        }

        for (my $i=0; $i< @{$newModel->{modelcompounds}}; $i++){
            my @msc1 = split /\//, $newModel->{modelcompounds}->[$i]->{compound_ref};
            my @msc = split /_/, $msc1[-1];

            if ($msc[0] eq 'cpd00000'){

                $eachModelCpds->{ $params->{output_model}}++;
              }
            else {

                $eachModelCpds->{ $params->{output_model}}++;
                $eachModelMSCpds->{ $params->{output_model}}++;
              }

        }

    };
    if ($@) {
           die "Error loading object from the workspace:\n".$@;
    }

    my $counterHash;
    foreach my $r (@newModelArr){

      if (exists $eachTemplateHash->{$r}){

        $counterHash->{$eachTemplateHash->{$r}->[1]}++;
      }
      else{

        $counterHash->{"ModelSEED"}++;

      }

    }

  my $htmlLink1 = "/kb/module/work/tmp/modelViz.html";
  open my $mData, ">", $htmlLink1  or die "Couldn't open modelViz file $!\n";
  print $mData qq{<html>
  <head>
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">
      google.charts.load("current", {packages:["corechart", "bar"]});
      google.charts.setOnLoadCallback(drawChart);
      function drawChart() {
        var data = google.visualization.arrayToDataTable([
          ['Organism_Name', 'Number of Reactions'],
          ['Aspergillus_terreus_NIH2624', $counterHash->{'Aspergillus_terreus_NIH2624'}],
          ['Candida_tropicali_MYA-3404', $counterHash->{'Candida_tropicali_MYA-3404'}],
          ['Candida_glabrata_ASM254', $counterHash->{'Candida_glabrata_ASM254'}],
          ['Kluyveromyces_lactis_NRRL', $counterHash->{'Kluyveromyces_lactis_NRRL'}],
          ['Saccharomyces_cerevisiae_5288c',$counterHash->{'Saccharomyces_cerevisiae_5288c'}],
          ['Neurospora_crassa_OR74A', $counterHash->{'Neurospora_crassa_OR74A'}]

        ]);

        var options = {
          title: 'Published Model Intergreation Statistics',
          is3D: true,
        };

        var chart = new google.visualization.PieChart(document.getElementById('piechart_3d'));
        chart.draw(data, options);
      }

       google.charts.setOnLoadCallback(drawBarColorsR);

       google.charts.setOnLoadCallback(drawBarColorsC);

      function drawBarColorsR() {
        var data = google.visualization.arrayToDataTable([
          ['Organism_Name', 'Total Number of Reactions', 'Number of ModelSEED Reactions'],
          ['Aspergillus_terreus_NIH2624', $eachModelRxns->{'Aspergillus_terreus_NIH2624'}, $eachModelMSRxns->{'Aspergillus_terreus_NIH2624'}],
          ['Candida_tropicali_MYA-3404', $eachModelRxns->{'Candida_tropicali_MYA-3404'}, $eachModelMSRxns->{'Candida_tropicali_MYA-3404'}],
          ['Candida_glabrata_ASM254', $eachModelRxns->{'Candida_glabrata_ASM254'}, $eachModelMSRxns->{'Candida_glabrata_ASM254'} ],
          ['Saccharomyces_cerevisiae_5288c',$eachModelRxns->{'Saccharomyces_cerevisiae_5288c'}, $eachModelMSRxns->{'Saccharomyces_cerevisiae_5288c'}],
          ['Neurospora_crassa_OR74A', $eachModelRxns->{'Neurospora_crassa_OR74A'}, $eachModelMSRxns->{'Neurospora_crassa_OR74A'}],
          ['Kluyveromyces_lactis_NRRL', $eachModelRxns->{'Kluyveromyces_lactis_NRRL'}, $eachModelMSRxns->{'Kluyveromyces_lactis_NRRL'}],
          ['$params->{output_model}', $eachModelRxns->{$params->{output_model}}, $eachModelMSRxns->{$params->{output_model}}]

        ]);

        var options = {
        title: 'Reaction intergration statiscs',
        chartArea: {width: '50%'},
        colors: ['#62f442', '#41a0f4'],
        hAxis: {
          title: 'Reactions: published models Vs user created model',
          minValue: 0,
           textStyle: {
            bold: true,
            fontSize: 14,
            color: '#ffab91'
          },
          titleTextStyle: {
            bold: true,
            fontSize: 12,
            color: '#4d4d4d'
          }
        },
        vAxis: {
          title: 'Organisms'
        }
      };
      var chart = new google.visualization.BarChart(document.getElementById('chart_rxns'));
      chart.draw(data, options);
    }

      function drawBarColorsC() {
        var data = google.visualization.arrayToDataTable([
          ['Organism_Name', 'Total Number of Compounds', 'Number of ModelSEED Compounds'],
          ['Aspergillus_terreus_NIH2624', $eachModelCpds->{'Aspergillus_terreus_NIH2624'}, $eachModelMSCpds->{'Aspergillus_terreus_NIH2624'}],
          ['Candida_tropicali_MYA-3404', $eachModelCpds->{'Candida_tropicali_MYA-3404'}, $eachModelMSCpds->{'Candida_tropicali_MYA-3404'}],
          ['Candida_glabrata_ASM254', $eachModelCpds->{'Candida_glabrata_ASM254'}, $eachModelMSCpds->{'Candida_glabrata_ASM254'} ],
          ['Saccharomyces_cerevisiae_5288c',$eachModelCpds->{'Saccharomyces_cerevisiae_5288c'}, $eachModelMSCpds->{'Saccharomyces_cerevisiae_5288c'}],
          ['Neurospora_crassa_OR74A', $eachModelCpds->{'Neurospora_crassa_OR74A'}, $eachModelMSCpds->{'Neurospora_crassa_OR74A'}],
          ['Kluyveromyces_lactis_NRRL', $eachModelCpds->{'Kluyveromyces_lactis_NRRL'}, $eachModelMSCpds->{'Kluyveromyces_lactis_NRRL'}],
          ['$params->{output_model}', $eachModelCpds->{$params->{output_model}}, $eachModelMSCpds->{$params->{output_model}}]

        ]);

        var options = {
        title: 'Compound intergration statiscs',
        chartArea: {width: '50%'},
        colors: ['#f44185', '#41a0f4'],
        hAxis: {
          title: 'Compounds: published models Vs user created model',
          minValue: 0,
           textStyle: {
            bold: true,
            fontSize: 14,
            color: '#ffab91'
          },
          titleTextStyle: {
            bold: true,
            fontSize: 12,
            color: '#4d4d4d'
          }
        },
        vAxis: {
          title: 'Organisms'
        }
      };
      var chart = new google.visualization.BarChart(document.getElementById('chart_cpds'));
      chart.draw(data, options);
    }
    </script>
  </head>
  <body>
    <div id="piechart_3d" style="width: 1200px; height: 500px;"></div>
    <div id="chart_rxns" style="width: 1200px; height: 500px;"></div>
    <div id="chart_cpds" style="width: 1200px; height: 500px;"></div>
  </body>
</html>
                  };
  close $mData;

   my $htmlLinkHash1 = {
        path => $htmlLink1,
        name => 'Integrated Published Model Statistics',
        description => 'Integrated Published Model Statistics PieChart'
    };


=head
 my $htmlLink2 = "/kb/module/work/tmp/modelViz2.html";
  open my $cData, ">", $htmlLink2  or die "Couldn't open modelViz2 file $!\n";
  print $cData qq{<html>
  <head>
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">
      google.charts.load("current", {packages:["corechart","bar"]});
      google.charts.setOnLoadCallback(drawBarColorsR);

       google.charts.setOnLoadCallback(drawBarColorsC);

      function drawBarColorsR() {
        var data = google.visualization.arrayToDataTable([
          ['Organism_Name', 'Total Number of Reactions', 'Number of ModelSEED Reactions'],
          ['Aspergillus_terreus_NIH2624', $eachModelRxns->{'Aspergillus_terreus_NIH2624'}, $eachModelMSRxns->{'Aspergillus_terreus_NIH2624'}],
          ['Candida_tropicali_MYA-3404', $eachModelRxns->{'Candida_tropicali_MYA-3404'}, $eachModelMSRxns->{'Candida_tropicali_MYA-3404'}],
          ['Candida_glabrata_ASM254', $eachModelRxns->{'Candida_glabrata_ASM254'}, $eachModelMSRxns->{'Candida_glabrata_ASM254'} ],
          ['Saccharomyces_cerevisiae_5288c',$eachModelRxns->{'Saccharomyces_cerevisiae_5288c'}, $eachModelMSRxns->{'Saccharomyces_cerevisiae_5288c'}],
          ['Neurospora_crassa_OR74A', $eachModelRxns->{'Neurospora_crassa_OR74A'}, $eachModelMSRxns->{'Neurospora_crassa_OR74A'}],
          ['$params->{output_model}', $eachModelRxns->{$params->{output_model}}, $eachModelMSRxns->{$params->{output_model}}]

        ]);

        var options = {
        title: 'Reaction intergration statiscs',
        chartArea: {width: '50%'},
        colors: ['#62f442', '#41a0f4'],
        hAxis: {
          title: 'Reactions: published models Vs user created model',
          minValue: 0,
           textStyle: {
            bold: true,
            fontSize: 14,
            color: '#ffab91'
          },
          titleTextStyle: {
            bold: true,
            fontSize: 12,
            color: '#4d4d4d'
          }
        },
        vAxis: {
          title: 'Organisms'
        }
      };
      var chart = new google.visualization.BarChart(document.getElementById('chart_rxns'));
      chart.draw(data, options);
    }

      function drawBarColorsC() {
        var data = google.visualization.arrayToDataTable([
          ['Organism_Name', 'Total Number of Compounds', 'Number of ModelSEED Compounds'],
          ['Aspergillus_terreus_NIH2624', $eachModelCpds->{'Aspergillus_terreus_NIH2624'}, $eachModelMSCpds->{'Aspergillus_terreus_NIH2624'}],
          ['Candida_tropicali_MYA-3404', $eachModelCpds->{'Candida_tropicali_MYA-3404'}, $eachModelMSCpds->{'Candida_tropicali_MYA-3404'}],
          ['Candida_glabrata_ASM254', $eachModelCpds->{'Candida_glabrata_ASM254'}, $eachModelMSCpds->{'Candida_glabrata_ASM254'} ],
          ['Saccharomyces_cerevisiae_5288c',$eachModelCpds->{'Saccharomyces_cerevisiae_5288c'}, $eachModelMSCpds->{'Saccharomyces_cerevisiae_5288c'}],
          ['Neurospora_crassa_OR74A', $eachModelCpds->{'Neurospora_crassa_OR74A'}, $eachModelMSCpds->{'Neurospora_crassa_OR74A'}],
          ['$params->{output_model}', $eachModelCpds->{$params->{output_model}}, $eachModelMSCpds->{$params->{output_model}}]

        ]);

        var options = {
        title: 'Compound intergration statiscs',
        chartArea: {width: '50%'},
        colors: ['#f44185', '#41a0f4'],
        hAxis: {
          title: 'Compounds: published models Vs user created model',
          minValue: 0,
           textStyle: {
            bold: true,
            fontSize: 14,
            color: '#ffab91'
          },
          titleTextStyle: {
            bold: true,
            fontSize: 12,
            color: '#4d4d4d'
          }
        },
        vAxis: {
          title: 'Organisms'
        }
      };
      var chart = new google.visualization.BarChart(document.getElementById('chart_cpds'));
      chart.draw(data, options);
    }
    </script>
  </head>
  <body>
    <div id="chart_rxns" style="width: 1200px; height: 500px;"></div>
    <div id="chart_cpds" style="width: 1200px; height: 500px;"></div>
  </body>
</html>
                  };
  close $cData;

   my $htmlLinkHash2 = {
        path => $htmlLink2,
        name => 'Integrated Reaction and Compounds Statistics',
        description => 'Integrated Published Model Statistics Bar-charts'
    };

=cut

    print &Dumper ($counterHash);

    my $stat_string1= "Fungal model was built based based on proteome comparison $protCompId and produced the model $params->{output_model}\n The intergreation of published models are as follows\n ";
    #my $stat_string2 = "\nAspergillus_terreus\t$counterHash->{'Aspergillus_terreus_NIH2624'}\nCandida_tropicali_MYA-3404\t$counterHash->{'Candida_tropicali_MYA-3404'}\nCandida_glabrata_ASM254\t$counterHash->{'Candida_glabrata_ASM254'}\nSaccharomyces_cerevisiae_5288c\t$counterHash->{'Saccharomyces_cerevisiae_5288c'}\nNeurospora_crassa_OR74A\t$counterHash->{'Neurospora_crassa_OR74A'}";

    my $reporter_string = $stat_string1;
    my $uid = UUID::Random::generate;
    my $report_context = {
      message => $reporter_string,
      objects_created => [],
      workspace_name => $params->{workspace},
      direct_html_link_index => 0,
      warnings => [],
      html_links => [$htmlLinkHash1],
      file_links =>[],
      report_object_name => "Report"."modelpropagation"."-".UUID::Random::generate
    };


    my $report_response;

    eval {
      $report_response = $reportHandle->create_extended_report($report_context);
    };
    if ($@){
      print "Exception message: " . $@->{"message"} . "\n";
      print "JSONRPC code: " . $@->{"code"} . "\n";
      print "Method: " . $@->{"method_name"} . "\n";
      print "Client-side exception:\n";
      print $@;
      print "Server-side exception:\n";
      print $@->{"data"};
      die $@;
    }

    print "Report is generated: name and the ref as follows\n";
    print &Dumper ($report_response);
    my $report_out = {
      report_name => $report_response->{name},
      report_ref => $report_response->{ref}
    };
    return $report_out;

    #END build_fungal_model
    my @_bad_returns;
    (ref($output) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"output\" (value was \"$output\")");
    if (@_bad_returns) {
    	my $msg = "Invalid returns passed to build_fungal_model:\n" . join("", map { "\t$_\n" } @_bad_returns);
    	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
    							       method_name => 'build_fungal_model');
    }
    return($output);
}




=head2 build_fungal_template

  $output = $obj->build_fungal_template($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_fungalmodeling.fungalReferenceModelBuildInput
$output is a kb_fungalmodeling.fungalReferenceModelBuildOutput
fungalReferenceModelBuildInput is a reference to a hash where the following keys are defined:
	workspace has a value which is a string
	reference_genome has a value which is a string
	reference_model has a value which is a string
	genome_ws has a value which is a string
	model_ws has a value which is a string
fungalReferenceModelBuildOutput is a reference to a hash where the following keys are defined:
	master_template_model_ref has a value which is a string
	master_template_genome_ref has a value which is a string

</pre>

=end html

=begin text

$params is a kb_fungalmodeling.fungalReferenceModelBuildInput
$output is a kb_fungalmodeling.fungalReferenceModelBuildOutput
fungalReferenceModelBuildInput is a reference to a hash where the following keys are defined:
	workspace has a value which is a string
	reference_genome has a value which is a string
	reference_model has a value which is a string
	genome_ws has a value which is a string
	model_ws has a value which is a string
fungalReferenceModelBuildOutput is a reference to a hash where the following keys are defined:
	master_template_model_ref has a value which is a string
	master_template_genome_ref has a value which is a string


=end text



=item Description



=back

=cut

sub build_fungal_template
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to build_fungal_template:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'build_fungal_template');
    }

    my $ctx = $kb_fungalmodeling::kb_fungalmodelingServer::CallContext;
    my($output);
    #BEGIN build_fungal_template

    print &Dumper ($params);
    my $token=$ctx->token;
    my $provenance=$ctx->provenance;
    my $wshandle = Workspace::WorkspaceClient->new($self->{'workspace-url'},token=>$token);

    my $fbaO = new fba_tools::fba_toolsClient( $self->{'callbackURL'},
                                                            ( 'service_version' => 'release',
                                                              'async_version' => 'release',
                                                            )
                                                           );

    print &Dumper ($wshandle);
    my $ws = $params->{workspace};


my $model_list = ['25992/65', '25992/70' ,'25992/60', '25992/54', '25992/26', '25992/99'];

=head
#published models considered for the template

      iJL1454 => ['iJL1454', 'Aspergillus terreus iJL1454'],
      iNX804  => ['iNX804','Candida_glabrata_ASM254'],
      iCT646 => ['iCT646','Candida_tropicali_MYA-3404'],
      iOD907 => ['iOD907','GCF_000002515.2'],
      iJDZ836 => ['iJDZ836','Neurospora_crassa_OR74A'],
=cut

=head
    my $fungal_temp_community_model = $fbaO->merge_metabolic_models_into_community_model({
       fbamodel_id_list => $model_list,
       fbamodel_output_id => "FungalTemplate",
       workspace => $params->{workspace},
       mixed_bag_model => 1

    });

    print &Dumper ($fungal_temp_community_model);
    die;
=cut
    my $crassaModel =  '25857/11/2';
    my $start_genome_id  = '25857/2/3'; # 'Neurospora_crassa_OR74A',
    my $start_model_name = '25857/11/2'; # iJDZ836";
    my $master_temp = '26394/2/1';#25857/25/1';
    my $first_model;


    my $masterBio;
    eval {
       #$masterBio = $wshandle->get_objects([{ref=>'25992/78'}])->[0]{data}{biomasses}->[0]{biomasscompounds};
       $masterBio = $wshandle->get_objects([{ref=>'25992/78'}])->[0]{data};

    };
    if ($@) {
       die "Error loading object from the workspace:\n".$@;
    }

print &Dumper ($masterBio);
die;


   my $biomass_cpd_remove = {
       biomass_id => 'bio1',
       biomass_compound_id => '',
       biomass_coefficient => 0
    };

    my $tempbiomassArr;
    for (my $i=0; $i< @{$masterBio}; $i++){

      print "$masterBio->[$i]->{modelcompound_ref}\n";

      my @cpdid = split /\//, $masterBio->[$i]->{modelcompound_ref};
      if ($cpdid[-1] =~ /cpd/){
        print "$cpdid[-1]\n";

      }
      else{
        $biomass_cpd_remove = {
           biomass_id => 'bio1',
           biomass_compound_id => $cpdid[-1],
           biomass_coefficient => 0
        };
        push (@{$tempbiomassArr},$biomass_cpd_remove );
      }

    }

    print &Dumper ($tempbiomassArr);

    my $edited_model = $fbaO->edit_metabolic_model({

        fbamodel_id => 'FungalTemplate',
        fbamodel_output_id => "master_fungal_template",
        workspace =>  $params->{workspace},
        compounds_to_add => [],
        compounds_to_change => [],
        biomasses_to_add => [],
        biomass_compounds_to_change => $tempbiomassArr,
        reactions_to_remove => [],
        reactions_to_change => [],
        reactions_to_add => [],
        edit_compound_stoichiometry => []


    });

    eval {
       $masterBio = $wshandle->get_objects([{ref=>$edited_model->{new_fbamodel_ref}}])->[0]{data};
       #$masterBio = $wshandle->get_objects([{ref=>$crassaModel}])->[0]{data}{biomasses};
    };
    if ($@) {
       die "Error loading object from the workspace:\n".$@;
    }

    #print &Dumper ($masterBio);

    $masterBio->{type} = "FungalGenomeScale";

    print "\n************ $masterBio->{type}\n";


    my $obj_info_list = undef;
    eval {
        $obj_info_list = $wshandle->save_objects({
            'workspace'=>$params->{workspace},
            'objects'=>[{
                'type'=>'KBaseFBA.FBAModel',
                'data'=> $masterBio,
                'name'=>'master_fungal_template',
                'provenance'=>$provenance
            }]
        });
    };
    if ($@) {
        die "Error saving modified genome object to workspace:\n".$@;
    }


    print &Dumper ($obj_info_list);

    die;


    #END build_fungal_template
    my @_bad_returns;
    (ref($output) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"output\" (value was \"$output\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to build_fungal_template:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'build_fungal_template');
    }
    return($output);
}




=head2 status

  $return = $obj->status()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module status. This is a structure including Semantic Versioning number, state and git info.

=back

=cut

sub status {
    my($return);
    #BEGIN_STATUS
    $return = {"state" => "OK", "message" => "", "version" => $VERSION,
               "git_url" => $GIT_URL, "git_commit_hash" => $GIT_COMMIT_HASH};
    #END_STATUS
    return($return);
}

=head1 TYPES



=head2 fungalmodelbuiltInput

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace has a value which is a string
genome_ref has a value which is a string
template_model has a value which is a string
gapfill_model has a value which is an int
translation_policy has a value which is a string
output_model has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace has a value which is a string
genome_ref has a value which is a string
template_model has a value which is a string
gapfill_model has a value which is an int
translation_policy has a value which is a string
output_model has a value which is a string


=end text

=back



=head2 fungalmodelbuiltOutput

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string


=end text

=back



=head2 fungalReferenceModelBuildInput

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace has a value which is a string
reference_genome has a value which is a string
reference_model has a value which is a string
genome_ws has a value which is a string
model_ws has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace has a value which is a string
reference_genome has a value which is a string
reference_model has a value which is a string
genome_ws has a value which is a string
model_ws has a value which is a string


=end text

=back



=head2 fungalReferenceModelBuildOutput

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
master_template_model_ref has a value which is a string
master_template_genome_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
master_template_model_ref has a value which is a string
master_template_genome_ref has a value which is a string


=end text

=back



=cut

1;
