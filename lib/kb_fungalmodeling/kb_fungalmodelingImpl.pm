package kb_fungalmodeling::kb_fungalmodelingImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org
our $VERSION = '0.0.1';
our $GIT_URL = 'https://github.com/janakagithub/kb_fungalmodeling.git';
our $GIT_COMMIT_HASH = 'c9c4015fce083ace29122468e7a63f8d0b2c7a28';

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





sub generate_pie_chart
{
    my($uMcounterHashGPR, $templateId) = @_;

    #Setting the pieChart on rxns
    my $piechartStringHeader = "['Organism_Name', 'Number of Reactions']";
    my $piechartString =[];
    push (@{$piechartString}, $piechartStringHeader);

    foreach my $k (keys $templateId){

        if ($k eq 'default_temp'){
            next;
        }
        else{
            my $tempStr = "['$templateId->{$k}->[1]',$uMcounterHashGPR->{$templateId->{$k}->[1]}]";
            push (@{$piechartString}, $tempStr);
        }

    }

    my $pieChartRxn = join(',' , @{$piechartString});

    return $pieChartRxn
}


sub generate_rxn_barchart
{
    my($eachModelRxns, $eachModelMSRxns, $eachModelGARxnsCount,$templateId, $uM, $uMgprRxnCount) = @_;
    my $barChartRxnsHeader = "['Organism_Name', 'Total Number of Reactions', 'GPR Reactions', 'Number of ModelSEED Reactions']";
    my $barChartRxnString = [];

    push (@{$barChartRxnString}, $barChartRxnsHeader );
    foreach my $k (keys $templateId){

        if ($k eq 'default_temp'){
            next;
        }
        else{
            my $tempStr = "['$templateId->{$k}->[1]',$eachModelRxns->{$templateId->{$k}->[1]},$eachModelGARxnsCount->{$templateId->{$k}->[1]},$eachModelMSRxns->{$templateId->{$k}->[1]} ]";
            push (@{$barChartRxnString}, $tempStr);
        }

    }
    my $uMstring = "['$uM', $eachModelRxns->{$uM}, $uMgprRxnCount, $eachModelMSRxns->{$uM}]";
    push (@{$barChartRxnString}, $uMstring);
    my $barString = join(',' , @{$barChartRxnString});
    return $barString;

}

sub generate_cpd_barchart
{
     my($eachModelCpds, $eachModelMSCpds,$templateId,$uM) = @_;
    my $barChartCpdsHeader = "['Organism_Name', 'Total Number of Compounds', 'Number of ModelSEED Compounds']";
    my $barChartCpdString = [];

    push (@{$barChartCpdString}, $barChartCpdsHeader );
    foreach my $k (keys $templateId){

        if ($k eq 'default_temp'){
            next;
        }
        else{
            my $tempStr = "['$templateId->{$k}->[1]',$eachModelCpds->{$templateId->{$k}->[1]},$eachModelMSCpds->{$templateId->{$k}->[1]} ]";
            push (@{$barChartCpdString}, $tempStr);
        }

    }
    my $uMstring = "['$uM', $eachModelCpds->{$uM}, $eachModelMSCpds->{$uM}]";
    push (@{$barChartCpdString}, $uMstring);

    my $barString = join(',' , @{$barChartCpdString});
    return $barString;

}


sub make_report_page
{

  my($pieChartRxn,$barChartRxn,$barChartCpds) = @_;
  my $htmlLink1 = "/kb/module/work/tmp/modelViz.html";
  open my $mData, ">", $htmlLink1  or die "Couldn't open modelViz file $!\n";
  print $mData qq{<html>
  <head>
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">
      google.charts.load("current", {packages:["corechart", "bar"]});
      google.charts.setOnLoadCallback(drawChart);
      function drawChart() {
        var data = google.visualization.arrayToDataTable([$pieChartRxn]);

        var options = {
          title: 'Published Model Integration Statistics',
          is3D: true,
        };

        var chart = new google.visualization.PieChart(document.getElementById('piechart_3d'));
        chart.draw(data, options);
      }

       google.charts.setOnLoadCallback(drawBarColorsR);
       google.charts.setOnLoadCallback(drawBarColorsC);

      function drawBarColorsR() {
        var data = google.visualization.arrayToDataTable([$barChartRxn]);

        var options = {
        title: 'Reaction integration statiscs',
        chartArea: {width: '50%'},
        colors: ['#62f442', '#41a0f4', '#ffab91'],
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
        var data = google.visualization.arrayToDataTable([$barChartCpds]);

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

   return $htmlLinkHash1;
}

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
	media_ref has a value which is a string
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
	media_ref has a value which is a string
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
    my $token=$ctx->token;
    my $provenance=$ctx->provenance;
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

    my $wshandle= Workspace::WorkspaceClient->new($self->{'workspace-url'},token=>$token);

    my $template_ws = 'jplfaria:narrative_1510597445008';# 'janakakbase:narrative_1513399583946'; # template workspaces
    my $template_genome_ref = 'FungalTemplate.genome';
    my $template_model_ref = 'master_fungal_template_fix';
    my $ws_name = $params->{workspace};
    my $protCompId = 'proteinComp'.$params->{genome_ref};
    my $tmpGenome;
    my $tmpModel;

    my $templateId = {
      default_temp => [$template_model_ref, $template_genome_ref,$template_genome_ref],
      iJL1454 => ['iJL1454_KBase', 'Aspergillus_terreus_NIH2624','GCF_000149615.1'],
      iNX804  => ['iNX804_KBase','Candida_glabrata_ASM254','GCF_000002545.3'],
      iCT646 => ['iCT646_KBase','Candida_tropicali_MYA-3404','GCF_000006335.2'],
      iOD907 => ['iOD907_KBase','Kluyveromyces_lactis_NRRL','GCF_000146045.2'],
      iJDZ836 => ['iJDZ836_KBase','Neurospora_crassa_OR74A','GCF_000182925.2'],
      iLC915 =>  ['iLC915_KBase', 'Komagataella_phaffii_GS115','GCF_000027005.1'],
      iRL766 => ['iRL766_KBase', 'Eremothecium_gossypii_ATCC_10895','GCF_000091025.4'],
      iAL1006 => ['iAL1006_KBase', 'Penicillium_rubens_Wisconsin','GCF_000226395.1'],
      iSS884 =>  ['iSS884_KBase', 'Scheffersomyces_stipitis_CBS','GCF_000209165.1'],
      iNL895 =>  ['iNL895_KBase', 'Yarrowia_lipolytica_CLIB122','GCF_000002525.2'],
      iWV1213 => ['iWV1213_KBase', 'Mucor_circinelloides_CBS277','Mucor_circinelloides_CBS277.49_v2.0'],
      iWV1314 => ['iWV1314_KBase', 'Aspergillus_oryzae_RIB40','GCF_000184455.2'],
      iMM904 => ['iMM904_KBase','Saccharomyces_cerevisiae_5288c','GCF_000146045.2']


    };

  my $eachTemplateHash;
  my $eachTemplateHashSplit;
  my $eachModelRxns;
  my $eachModelMSRxns;
  my $eachModelCpds;
  my $eachModelMSCpds;
  my $eachModelGARxns;
  my $eachModelGARCount;
  my $eachModelNonGARxns;
  my $eachModelBM;
  my $eachModelGARxnsCount;
  my @newModelArr;
  my $templateGenomeRefs;

# Generate stats
foreach my $k (keys $templateId){
    if ($k eq 'default_temp'){
      $templateGenomeRefs->{$k} = '25992/104';
      next;
    }
    else{
        eval {
            print "retrieving individual models from template $k\n";
            my $eachTemplate = $wshandle->get_objects([{workspace=>$template_ws,name=>$templateId->{$k}->[0]}] )->[0]{data};#{modelreactions};
            $templateGenomeRefs->{$k} = $eachTemplate->{genome_ref};
            $eachModelBM->{$templateId->{$k}->[1]} = $eachTemplate->{biomasses};
            for (my $i=0; $i< @{$eachTemplate->{modelreactions}}; $i++){

                my @rid = split /_/, $eachTemplate->{modelreactions}->[$i]->{id};
                $eachTemplateHashSplit->{ $templateId->{$k}->[1] }->{ $rid[0] }  = [$k, $templateId->{$k}->[1]];

                $eachTemplateHash->{ $templateId->{$k}->[1] }->{ $eachTemplate->{modelreactions}->[$i]->{id}}  = [$k, $templateId->{$k}->[1]];
                my @msr1 = split /\//, $eachTemplate->{modelreactions}->[$i]->{reaction_ref};
                my @msr = split /_/, $msr1[-1];

                if ($msr[0] eq 'rxn00000'){

                    $eachModelRxns->{$templateId->{$k}->[1]}++;

                    if ($eachTemplate->{modelreactions}->[$i]->{imported_gpr}){

                        $eachModelGARxns->{ $templateId->{$k}->[1] }->{$eachTemplate->{modelreactions}->[$i]->{id}} = [$eachTemplate->{modelreactions}->[$i]->{imported_gpr},$eachTemplate->{modelreactions}->[$i]->{name}];
                        $eachModelGARxnsCount->{$templateId->{$k}->[1]}++;
                    }
                    else{
                        push (@{$eachModelNonGARxns->{$templateId->{$k}->[1]}}, $eachTemplate->{modelreactions}->[$i]->{id});
                    }


                }
                else {

                      $eachModelRxns->{$templateId->{$k}->[1]}++;
                      $eachModelMSRxns->{$templateId->{$k}->[1]}++;

                      if ($eachTemplate->{modelreactions}->[$i]->{imported_gpr}){

                        $eachModelGARxns->{ $templateId->{$k}->[1] }->{$eachTemplate->{modelreactions}->[$i]->{id}} = [$eachTemplate->{modelreactions}->[$i]->{imported_gpr},$eachTemplate->{modelreactions}->[$i]->{name}];
                        $eachModelGARxnsCount->{$templateId->{$k}->[1]}++;

                      }
                      else{
                        push (@{$eachModelNonGARxns->{$templateId->{$k}->[1]}}, $eachTemplate->{modelreactions}->[$i]->{id});

                      }

                }

            } # for

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
  #print &Dumper ($eachModelRxns);
  #print &Dumper ($eachModelMSRxns);
  #print &Dumper ($eachModelCpds);
  #print &Dumper ($eachModelMSCpds);

    if (defined $params->{template_model}) {

      print "Template selected as $params->{template_model} \n";
      $tmpModel = $templateId->{$params->{template_model}}->[0];
      $tmpGenome = $templateId->{$params->{template_model}}->[2];
      #$tmpGenome = $templateGenomeRefs->{$params->{template_model}};
      print &Dumper ($templateId);

    }
    my $tran_policy;
    if (defined $params->{translation_policy}){
      $tran_policy = $params->{translation_policy};

    }

#=head # commenting the protein comparison
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
#=cut

    my $gpModelFromSource;
    if ($params->{gapfill_model} == 1){

        my $dr_model =$params->{genome_ref}."_draftModel";
        my $fba_modelProp =  $fbaO->propagate_model_to_new_genome({
            fbamodel_id => $tmpModel,
            fbamodel_workspace => $template_ws,
            proteincomparison_id => $protCompId, #'proteinCompPsean1',
            proteincomparison_workspace => $params->{workspace},
            fbamodel_output_id =>  $dr_model,
            workspace => $params->{workspace},
            keep_nogene_rxn => 0,
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
            #source_fbamodel_id => 'Fungi', #$tmpModel,
            #source_fbamodel_workspace => 'NewKBaseModelTemplates', #$template_ws,
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

#=head
        my $fba_modelProp =  $fbaO->propagate_model_to_new_genome({
            fbamodel_id => $tmpModel,
            fbamodel_workspace => $template_ws,
            proteincomparison_id => $protCompId, #'proteinCompPsean1'
            proteincomparison_workspace => $params->{workspace},
            fbamodel_output_id =>  $params->{output_model},
            workspace => $params->{workspace},
            keep_nogene_rxn => 0,
            #media_id =>
            #media_workspace =>
            minimum_target_flux => 0.1,
            translation_policy => $tran_policy

        });
#=cut
    }

     my $userModelRxns;
     my $newModel;
     #temperoraly fetching the object - remove after debugging
     #$params->{output_model} = 'prpogated_model_out_Psean1';
    eval {
        $newModel = $wshandle->get_objects([{workspace=>$params->{workspace}, name=>$params->{output_model}}])->[0]{data};# ->{modelreactions};

        for (my $i=0; $i< @{$newModel->{modelreactions}}; $i++){
            $userModelRxns->{$newModel->{modelreactions}->[$i]->{id}} =1;
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

    #########################\
    my $uMcounterHashGPR;
    my $uMcounterHashNOGPR;
    my $uMGPRrxns;
    my $uMNOGPRrxns;
    my $uMgprRxnCount =0;
    for (my $i=0; $i< @{$newModel->{modelreactions}}; $i++){

        my $gprCheck = $newModel->{modelreactions}->[$i];
        my @rid = split /_/, $gprCheck->{id};
        if (!@{$gprCheck->{modelReactionProteins} }){

            foreach my $k (keys $templateId){
                if ($k eq 'default_temp'){
                    next;
                }
                elsif (exists $eachTemplateHash->{ $templateId->{$k}->[1] }->{ $gprCheck->{id} }){
                    $uMcounterHashNOGPR->{$templateId->{$k}->[1]}++;
                    $uMNOGPRrxns->{$gprCheck->{id}} = 1;
                    #print "NO GPR $gprCheck->{id}\t $templateId->{$k}->[1] \n";
                }
                elsif (exists $eachTemplateHashSplit->{ $templateId->{$k}->[1] }->{ $rid[0] }){
                    $uMcounterHashNOGPR->{$templateId->{$k}->[1]}++;
                    $uMNOGPRrxns->{$gprCheck->{id}} = 1;
                    #print "NO GPR $gprCheck->{id}\t $templateId->{$k}->[1] \n";
                }
                else{
                  next;
                }

            }

        }

        else{

            my $gprFlag =0;
            my $eachGprArr = [];
            for (my $j=0; $j< @{$gprCheck->{modelReactionProteins}}; $j++){

                my $subGPR = $gprCheck->{modelReactionProteins}->[$j];

                for (my $n=0; $n< @{$subGPR->{modelReactionProteinSubunits}}; $n++){
                    my $eachGPR = $subGPR->{modelReactionProteinSubunits}->[$n];

                    if (!@{$eachGPR->{feature_refs} }){
                        next;

                    }
                    else{

                        $gprFlag =1;
                        $eachGprArr = $eachGPR->{feature_refs};
                        last;
                    }
                }
            }
            if ($gprFlag != 0){
                foreach my $k (keys $templateId){
                    if ($k eq 'default_temp'){
                    next;
                    }
                    elsif (exists $eachTemplateHash->{ $templateId->{$k}->[1] }->{ $gprCheck->{id} }){
                        $uMcounterHashGPR->{$templateId->{$k}->[1]}++;
                        $uMGPRrxns->{$gprCheck->{id}} = $eachGprArr;
                        $uMgprRxnCount++;
                        #print "GPR $gprCheck->{id}\t $templateId->{$k}->[1] \n";
                    }
                    elsif (exists $eachTemplateHashSplit->{ $templateId->{$k}->[1] }->{ $rid[0] }){
                        $uMcounterHashGPR->{$templateId->{$k}->[1]}++;
                        $uMGPRrxns->{$gprCheck->{id}} = $eachGprArr;
                        $uMgprRxnCount++
                        #print "GPR $gprCheck->{id}\t $templateId->{$k}->[1] \n";
                    }
                    else{
                      next;
                    }
                }
            }
            else{
                foreach my $k (keys $templateId){
                    if ($k eq 'default_temp'){
                    next;
                    }
                    elsif (exists $eachTemplateHash->{ $templateId->{$k}->[1] }->{ $gprCheck->{id} }){
                        $uMcounterHashNOGPR->{$templateId->{$k}->[1]}++;
                        $uMNOGPRrxns->{$gprCheck->{id}} = 1;
                        #print "NO GPR $gprCheck->{id}\t $templateId->{$k}->[1] \n";
                    }
                    elsif (exists $eachTemplateHashSplit->{ $templateId->{$k}->[1] }->{ $rid[0] }){
                        $uMcounterHashNOGPR->{$templateId->{$k}->[1]}++;
                        $uMNOGPRrxns->{$gprCheck->{id}} = 1;
                        #print "NO GPR $gprCheck->{id}\t $templateId->{$k}->[1] \n";
                    }
                    else{
                      next;
                    }
                }
            }
        }

    }

    my @sortedUmCounterHash;
    @sortedUmCounterHash = sort {$uMcounterHashGPR->{$a} <=> $uMcounterHashGPR->{$b}} keys (%$uMcounterHashGPR);
    #print &Dumper ($uMcounterHashGPR);
    #print &Dumper ($uMcounterHashNOGPR);
    #print &Dumper (\@sortedUmCounterHash);
    #print &Dumper ($uMNOGPRrxns);

    ########################### compute the closeset template model and prep a list of non-gpr reactions to remove
    my $removeRxnsArr = [];
    if (defined $uMNOGPRrxns){
        foreach my $k (keys $uMNOGPRrxns){
            my @rid = split /_/, $k;
            if (@sortedUmCounterHash && exists $eachTemplateHash->{ $sortedUmCounterHash[-1] }->{ $k }){
                next;
            }
            elsif (@sortedUmCounterHash && exists $eachTemplateHashSplit->{ $sortedUmCounterHash[-1] }->{ $rid[0] }){
                next;
            }
            else {

                push (@{$removeRxnsArr}, $k);
            }
        }
    }
    foreach my $k (@{$eachModelNonGARxns->{$sortedUmCounterHash[-1] } }){
        #print "$k\t". scalar(@{$eachModelNonGARxns->{$sortedUmCounterHash[-1] } })."\n";
    }
    #Removing reactions from compartments if the reaction count is less than 15

    my $compartmentRxnHash = {};
    for (my $i=0; $i< @{$newModel->{modelreactions}}; $i++){

      my @compt = split /_/, $newModel->{modelreactions}->[$i]->{id};

        #print "$compt[-1]\t$newModel->{modelreactions}->[$i]->{id}\n";
      #if ($compt[-1] eq 'n0' || $compt[-1] eq 'g0' || $compt[-1] eq 'v0' || $compt[-1] eq 'p0' || $compt[-1] eq 'z0' ||$compt[-1] eq 'r0' ){
        push (@{$compartmentRxnHash->{$compt[-1]}}, $newModel->{modelreactions}->[$i]->{id});
      #}

    }

    my $compt_reaction_count = 15;
    my $compt_reaction_hash;
    my $compt_reaction_hash_mixed;
    my $compt_Z_reaction_arr = [];

    foreach my $k (keys $compartmentRxnHash){
      if ($k eq 'e0' || $k eq 'c0' || $k eq 'm0' || $k eq 'x0'){

        for (my $i=0; $i< @{$compartmentRxnHash->{$k}}; $i++){

          $compt_reaction_hash->{$compartmentRxnHash->{$k}->[$i]} = $k;
          my @tmprxn = split /_/, $compartmentRxnHash->{$k}->[$i];
          $compt_reaction_hash_mixed->{$tmprxn[0]} = $k;
        }

        next;
      }
      elsif ( $k eq 'z0'){
        for (my $i=0; $i< @{$compartmentRxnHash->{$k}}; $i++){
          push (@{$compt_Z_reaction_arr}, $compartmentRxnHash->{$k}->[$i]);
          push (@{$removeRxnsArr}, $compartmentRxnHash->{$k}->[$i]);
        }

      }
      else{

          if (@{$compartmentRxnHash->{$k}} < $compt_reaction_count){

            for (my $i=0; $i< @{$compartmentRxnHash->{$k}}; $i++){
              push (@{$removeRxnsArr}, $compartmentRxnHash->{$k}->[$i]);
            }

          }
          else{
            for (my $i=0; $i< @{$compartmentRxnHash->{$k}}; $i++){
                $compt_reaction_hash->{$compartmentRxnHash->{$k}->[$i]} = $k;
                my @tmprxn = split /_/, $compartmentRxnHash->{$k}->[$i];
                $compt_reaction_hash_mixed->{$tmprxn[0]} = $k;
            }

          }
      }

    }

    print &Dumper ($uMGPRrxns);


    #Adding the Z compartment reactions as cytosol reactions
    #print &Dumper ($compt_Z_reaction_arr);
    my $addReactionArr = [];
    for (my $i=0; $i< @{$compt_Z_reaction_arr}; $i++){

        if ($compt_Z_reaction_arr->[$i] =~ /rxn/){
            my $comptFlag = 0;

                my @temprxn = split /_/, $compt_Z_reaction_arr->[$i];

                if (exists $compt_reaction_hash_mixed->{$temprxn[0]}){
                    print "$compt_Z_reaction_arr->[$i]\n";
                    $comptFlag = 1;

                }
                else{
                    my @gprarr = split /\//, $uMGPRrxns->{$compt_Z_reaction_arr->[$i]}->[0];
                    my $newRxntoAdd = {
                        add_reaction_id => $temprxn[0],
                        reaction_compartment_id => ["c0"],
                        add_reaction_direction => '<=>',
                        add_reaction_gpr => $gprarr[-1]

                    };
                    push (@{$addReactionArr}, $newRxntoAdd);

                }
        }

    }

    print &Dumper ($addReactionArr);


    ## Generating Biomass
    my $tempbiomassArr;
    my $biomass_cpd_remove = {
       biomass_id => '',
       biomass_compound_id => '',
       biomass_coefficient => 0
    };


    for (my $j=0; $j< @{$eachModelBM->{ $sortedUmCounterHash[-1]} }; $j++){
        my $currentBiomass = $eachModelBM->{$sortedUmCounterHash[-1]}->[$j];
        for (my $i=0; $i< @{$currentBiomass->{biomasscompounds}}; $i++){

            my @cpdid = split /\//, $currentBiomass->{biomasscompounds}->[$i]->{modelcompound_ref};
            if ($cpdid[-1] =~ /cpd/){
                print "$cpdid[-1]\n";

                 $biomass_cpd_remove = {
                   biomass_id => $currentBiomass->{id},
                   biomass_compound_id => $cpdid[-1],
                   biomass_coefficient => $currentBiomass->{biomasscompounds}->[$i]->{coefficient}
                };
                push (@{$tempbiomassArr},$biomass_cpd_remove );

            }
            else{
                $biomass_cpd_remove = {
                   biomass_id => $currentBiomass->{id},
                   biomass_compound_id => $cpdid[-1],
                   biomass_coefficient => $currentBiomass->{biomasscompounds}->[$i]->{coefficient}
                };
                push (@{$tempbiomassArr},$biomass_cpd_remove );
            }
        }
    }

    print &Dumper ($removeRxnsArr);

    my $edited_model = $fbaO->edit_metabolic_model({

        fbamodel_id => $params->{workspace}.'/'.$params->{output_model},
        fbamodel_output_id => $params->{output_model},
        workspace =>  $params->{workspace},
        compounds_to_add => [],
        compounds_to_change => [],
        biomasses_to_add => [],
        biomass_compounds_to_change => [],
        reactions_to_remove => join (',',@{$removeRxnsArr}),
        reactions_to_change => [],
        reactions_to_add => $addReactionArr,
        edit_compound_stoichiometry => []


    });


    my $pieChartRxn = generate_pie_chart ($uMcounterHashGPR, $templateId);
    my $barChartRxn = generate_rxn_barchart($eachModelRxns, $eachModelMSRxns, $eachModelGARxnsCount,$templateId,$params->{output_model},$uMgprRxnCount );
    my $barChartCpds = generate_cpd_barchart($eachModelCpds, $eachModelMSCpds,$templateId,$params->{output_model});
    my $htmlLink1 = make_report_page ($pieChartRxn,$barChartRxn,$barChartCpds);
    #print $pieChartRxn ."\n".$barChartRxn. "\n". $barChartCpds ."\n";


    my $stat_string1= "Fungal model was built based based on proteome comparison $protCompId and produced the model $params->{output_model}\n The intergreation of published models are as follows\n ";
    #my $stat_string2 = "\nAspergillus_terreus\tuMcounterHashGPR0->{'Aspergillus_terreus_NIH2624'}\nCandida_tropicali_MYA-3404\tuMcounterHashGPR0->{'Candida_tropicali_MYA-3404'}\nCandida_glabrata_ASM254\tuMcounterHashGPR0->{'Candida_glabrata_ASM254'}\nSaccharomyces_cerevisiae_5288c\tuMcounterHashGPR0->{'Saccharomyces_cerevisiae_5288c'}\nNeurospora_crassa_OR74A\tuMcounterHashGPR0->{'Neurospora_crassa_OR74A'}";
    my $reporter_string = $stat_string1;
    my $uid = UUID::Random::generate;
    my $report_context = {
      message => $reporter_string,
      objects_created => [],
      workspace_name => $params->{workspace},
      direct_html_link_index => 0,
      warnings => [],
      html_links => [$htmlLink1],
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
my $model_listExtended = ["25992/65/4", "25992/70/4", "25992/60/2", "25992/54/3", "25992/26/10", "25992/99/2", "25992/130/1", "25992/139/1",  "25992/155/1", "25992/158/1", "25992/178/1", "25992/175/1"];

#merging into community model is done in the narrative.
=head
    my $fungal_temp_community_model = $fbaO->merge_metabolic_models_into_community_model({
       fbamodel_id_list => $model_list,
       fbamodel_output_id => "FungalTemplate",
       workspace => $params->{workspace},
       mixed_bag_model => 1

    });

    print &Dumper ($fungal_temp_community_model);
    die;


#build a microbial model here with the new tmeplate to extract the biomass cpds and coefficients

 my $build_model = $fbaO->build_metabolic_model ({

            genome_id => 'GCF_000146045.2',
            genome_workspace => 'jplfaria:narrative_1510597445008',
            workspace => 'jplfaria:narrative_1510597445008',
            template_id => 'Fungi',
            template_workspace => 'NewKBaseModelTemplates',
            fbamodel_output_id => 'ModelWithFungalBiomass'

        });


 my $fungalBiomassModel;
     eval {
       $fungalBiomassModel = $wshandle->get_objects([{ref=>'25992/232'}])->[0]{data}{modelreactions};
       #$masterBio = $wshandle->get_objects([{ref=>$crassaModel}])->[0]{data}{biomasses};
    };
    if ($@) {
       die "Error loading object from the workspace:\n".$@;
    }

  my $removeRxnBiomassModel = [];
 for (my $i=0; $i < @{$fungalBiomassModel}; $i++){
    #my @rid = split /_/, $fungalBiomassModel->[$i]->{id};
    push (@{$removeRxnBiomassModel},  $fungalBiomassModel->[$i]->{id});

 }

 my $edited_model = $fbaO->edit_metabolic_model({

        fbamodel_id => 'ModelWithFungalBiomass', # 'GF_GMM'
        fbamodel_output_id => 'ModelWithFungalBiomassEmpty', #'test_modifiedmodel'
        workspace =>  $params->{workspace}, #'janakakbase:narrative_1518190880851',
        compounds_to_add => [],
        compounds_to_change => [],
        biomasses_to_add => [],
        biomass_compounds_to_change => [],
        #biomass_compounds_to_change => [],
        reactions_to_remove => join (',',@{$removeRxnBiomassModel}),
        reactions_to_change => [],
        reactions_to_add => [],
        edit_compound_stoichiometry => []


    });



=cut
=head

# 1.) Here remove the existing biomass of the template model/ relink the Fungal template/ and change the type
    my $masterBio;
     eval {
       $masterBio = $wshandle->get_objects([{ref=>'25992/78'}])->[0]{data};
       #$masterBio = $wshandle->get_objects([{ref=>$crassaModel}])->[0]{data}{biomasses};
    };
    if ($@) {
       die "Error loading object from the workspace:\n".$@;
    }

    #print &Dumper ($masterBio);

    $masterBio->{biomasses} = [];
    $masterBio->{template_ref} = 'NewKBaseModelTemplates/Fungi';
    $masterBio->{type} = 'FungalGenomeScale';

    print "\n************ $masterBio->{type}\n";


    my $obj_info_list = undef;
    eval {
        $obj_info_list = $wshandle->save_objects({
            'workspace'=>$params->{workspace},
            'objects'=>[{
                'type'=>'KBaseFBA.FBAModel',
                'data'=> $masterBio,
                'name'=>'FungalTemplateEmpty',
                'provenance'=>$provenance
            }]
        });
    };
    if ($@) {
        die "Error saving modified genome object to workspace:\n".$@;
    }

# 2.) Adding new biomass to the template should be done in the Narrative : Use edit model

=cut

 # 3.) Test for biomass and other info by fetching the object
my $templateModel;
    eval {
       #$templateModel = $wshandle->get_objects([{ref=>'25992/228'}])->[0]{data}{biomasses}->[0]{biomasscompounds}; #Object Name : FungalTemplateBaseBiomass
       $templateModel = $wshandle->get_objects([{ref=>'25992/78'}])->[0]{data}{biomasses};#->[0]{biomasscompounds}; # Object Name : FungalTemplate

    };
    if ($@) {
       die "Error loading object from the workspace:\n".$@;
    }


    print &Dumper ($templateModel);
    die;




    my $crassaModel =  '25857/11/2';
    my $start_genome_id  = '25857/2/3'; # 'Neurospora_crassa_OR74A',
    my $start_model_name = '25857/11/2'; # iJDZ836";
    my $master_temp = '26394/2/1';#25857/25/1';
    my $first_model;


   my $biomass_cpd_remove = {
       biomass_id => 'bio1',
       biomass_compound_id => '',
       biomass_coefficient => 0
    };

    my $tempbiomassArr;
    for (my $i=0; $i< @{$templateModel}; $i++){

      print "$templateModel->[$i]->{modelcompound_ref}\t $i\n";

      my @cpdid = split /\//, $templateModel->[$i]->{modelcompound_ref};
      if ($cpdid[-1] =~ /cpd/){
        #print "$cpdid[-1]\n";


        ############################################
        #Biomass is removed from the master template as it complicates the gapfilling, instead individual model template biomass is intergrated
        #This can revert back, when we have sufficent level of modelSEED reaction/compound integration, which is not the case yet

         $biomass_cpd_remove = {
           biomass_id => 'bio1',
           biomass_compound_id => $cpdid[-1],
           biomass_coefficient => 0
        };
        push (@{$tempbiomassArr},$biomass_cpd_remove );


        #############################################

      }
      else{

        print "*****$cpdid[-1]********\n";
        die;
        print &Dumper ($templateModel->[$i]);
                $biomass_cpd_remove = {
           biomass_id => 'bio1',
           biomass_compound_id => $cpdid[-1],
           biomass_coefficient => 0
        };
        push (@{$tempbiomassArr},$biomass_cpd_remove );
      }

    }

    #print &Dumper ($templateModel);
    #die;

    #print &Dumper ($tempbiomassArr);

  my $biomassModel;
    eval {
       $biomassModel = $wshandle->get_objects([{ref=>'25992/235'}])->[0]{data}{biomasses}->[0]{biomasscompounds}; #Object Name : ModelWithFungalBiomassEmpty
       #$templateModel = $wshandle->get_objects([{ref=>'25992/83'}])->[0]{data}{biomasses};#->[0]{biomasscompounds}; #Object Name : master_fungal_template

    };
    if ($@) {
       die "Error loading object from the workspace:\n".$@;
    }

    my $addbiomass = [];
    my $biomass_cpd_add = {
       biomass_id => 'bio1',
       biomass_compound_id => '',
       biomass_coefficient => 0
    };

    for (my $i=0; $i< @{$biomassModel}; $i++){

      print "$templateModel->[$i]->{modelcompound_ref}\t $i\n";

      my @cpdid = split /\//, $biomassModel->[$i]->{modelcompound_ref};
      if ($cpdid[-1] =~ /cpd/){
        #print "$cpdid[-1]\n";


        my $coefficient = $biomassModel->[$i]->{coefficient} + 0;
         $biomass_cpd_add = {
           biomass_id => 'bio1',
           biomass_compound_id => $cpdid[-1],
           biomass_coefficient => $coefficient
        };
        push (@{$addbiomass},$biomass_cpd_add );
      }
    }

=head

    open my $fbm, "<", "/kb/module/data/BiomassCompounds.tsv" or die "Couldn't open Biomass file $!\n";

    <$fbm>;
    my $biomass_cpd_add = {
       biomass_id => 'bio1',
       biomass_compound_id => '',
       biomass_coefficient => 0
    };

    my $counter =0;
    my $tempbiomassAdd;
    while (my $input = <$fbm>){
        chomp $input;
        my @bcpdid = split /\t/, $input;

        print "**$bcpdid[1]**\t###$bcpdid[2]####\n";

        $bcpdid[2] = $bcpdid[2] + 0;

         $biomass_cpd_add = {
           biomass_id => 'bio1',
           biomass_compound_id => $bcpdid[1].'_c0',
           biomass_coefficient => $bcpdid[2]
        };
        push (@{$tempbiomassAdd},$biomass_cpd_add);


    }
=cut

    my $newbiomass = {
          biomass_name => 'bio1',
          biomass_dna => 0.004,
          biomass_rna => 0.0063,
          biomass_protein => 0.45,
          biomass_cellwall => 0.407,
          biomass_lipid => 0.029,
          biomass_cofactor => 0.047,
          biomass_energy => 60

        };


    my $edited_model = $fbaO->edit_metabolic_model({

        fbamodel_id => 'FungalTemplateBaseBiomass', # 'GF_GMM'
        fbamodel_output_id => 'FungalTemplateBaseBiomassEmpty', #'test_modifiedmodel'
        workspace =>  $params->{workspace}, #'janakakbase:narrative_1518190880851',
        compounds_to_add => [],
        compounds_to_change => [],
        biomasses_to_add => [],
        biomass_compounds_to_change => $tempbiomassArr,
        #biomass_compounds_to_change => [],
        reactions_to_remove => [],
        reactions_to_change => [],
        reactions_to_add => [],
        edit_compound_stoichiometry => []


    });

    my $edited_model = $fbaO->edit_metabolic_model({

        fbamodel_id => 'FungalTemplateBaseBiomassEmpty', # 'GF_GMM'
        fbamodel_output_id => 'master_fungal_template', #'test_modifiedmodel'
        workspace =>  $params->{workspace}, #'janakakbase:narrative_1518190880851',
        compounds_to_add => [],
        compounds_to_change => [],
        biomasses_to_add => [],
        biomass_compounds_to_change => $addbiomass,
        #biomass_compounds_to_change => [],
        reactions_to_remove => [],
        reactions_to_change => [],
        reactions_to_add => [],
        edit_compound_stoichiometry => []


    });
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




=head2 build_model_stats

  $output = $obj->build_model_stats($params)

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

sub build_model_stats
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to build_model_stats:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'build_model_stats');
    }

    my $ctx = $kb_fungalmodeling::kb_fungalmodelingServer::CallContext;
    my($output);
    #BEGIN build_model_stats
    my $token=$ctx->token;
    my $provenance=$ctx->provenance;
    my $output = "stats_return";
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

    my $wshandle= Workspace::WorkspaceClient->new($self->{'workspace-url'},token=>$token);


    my $template_ws = 'jplfaria:narrative_1510597445008';# 'janakakbase:narrative_1513399583946'; # template workspaces
    my $template_genome_ref = 'FungalTemplate.genome';
    my $template_model_ref = 'master_fungal_template';
    my $ws_name = $params->{workspace};
    my $protCompId = 'proteinComp'.$params->{genome_ref};
    my $tmpGenome;
    my $tmpModel;


    my $templateId = {
      #default_temp => [$template_model_ref, $template_genome_ref],
      iJL1454 => ['iJL1454_KBase', 'Aspergillus_terreus_NIH2624','GCF_000149615.1'],
      iNX804  => ['iNX804_KBase','Candida_glabrata_ASM254','GCF_000002545.3'],
      iCT646 => ['iCT646_KBase','Candida_tropicali_MYA-3404','GCF_000006335.2'],
      iOD907 => ['iOD907_KBase','Kluyveromyces_lactis_NRRL','GCF_000146045.2'],
      iJDZ836 => ['iJDZ836_KBase','Neurospora_crassa_OR74A','GCF_000182925.2'],
      iLC915 =>  ['iLC915_KBase', 'Komagataella_phaffii_GS115','GCF_000027005.1'],
      iRL766 => ['iRL766_KBase', 'Eremothecium_gossypii_ATCC_10895','GCF_000091025.4'],
      iAL1006 => ['iAL1006_KBase', 'Penicillium_rubens_Wisconsin','GCF_000226395.1'],
      iSS884 =>  ['iSS884_KBase', 'Scheffersomyces_stipitis_CBS','GCF_000209165.1'],
      iNL895 =>  ['iNL895_KBase', 'Yarrowia_lipolytica_CLIB122','GCF_000002525.2'],
      iWV1213 => ['iWV1213_KBase', 'Mucor_circinelloides_CBS277','Mucor_circinelloides_CBS277.49_v2.0'],
      iWV1314 => ['iWV1314_KBase', 'Aspergillus_oryzae_RIB40','GCF_000184455.2'],
      iMM904 => ['iMM904_KBase','Saccharomyces_cerevisiae_5288c','GCF_000146045.2']

    };
=head
    foreach my $k (keys $templateId){

        my $dr_model =$k."_DM";
        my $fbaMid = $k.'_CM';
        my $prId = $k.'_P';

        print $dr_model ."\t". $fbaMid . "\t". $prId ."\n";


        my $fba_modelProp =  $fbaO->propagate_model_to_new_genome({
            fbamodel_id => "28321/40/1", #$fbaMid,
            fbamodel_workspace => 'janakakbase:narrative_1518296212814',
            proteincomparison_id => "28321/62/1", #$prId,
            proteincomparison_workspace => 'janakakbase:narrative_1518296212814',
            fbamodel_output_id =>  $dr_model,
            workspace => 'janakakbase:narrative_1518296212814',
            keep_nogene_rxn => 1,
            #media_id =>
            #media_workspace =>
            minimum_target_flux => 0.1,
            translation_policy => 'translate_only'
            #output_id =>  $dr_model
        });
        die;
        print &Dumper ($fba_modelProp);

    }
    die;
=cut
      my $eachTemplateHash;
      my $eachTemplateHashSplit;
      my $eachModelRxns;
      my $eachModelMSRxns;
      my $eachModelCpds;
      my $eachModelMSCpds;
      my $eachModelGARxns;
      my $eachModelGARCount;
      my $eachModelMSGARxns;
      my $eachModelMSGARCount;
      my $eachModelNonGARxns;
      my $eachModelBM;
      my $eachModelGARxnsCount;
      my @newModelArr;
      my $templateGenomeRefs;



    # Generate stats
    foreach my $k (keys $templateId){
        if ($k eq 'default_temp'){
          #$templateGenomeRefs->{$k} = '25992/104';
          next;
        }
        else{
            eval {
                print "retrieving individual models from template $k\n";
                my $eachTemplate = $wshandle->get_objects([{workspace=>$template_ws,name=>$templateId->{$k}->[0]}] )->[0]{data};#{modelreactions};
                $templateGenomeRefs->{$k} = $eachTemplate->{genome_ref};
                $eachModelBM->{$templateId->{$k}->[0]} = $eachTemplate->{biomasses};

                for (my $i=0; $i< @{$eachTemplate->{modelreactions}}; $i++) {
                    my @rid = split /_/, $eachTemplate->{modelreactions}->[$i]->{id};

                    my @msr1 = split /\//, $eachTemplate->{modelreactions}->[$i]->{reaction_ref};
                    my @msr = split /_/, $msr1[-1];

                    $eachTemplateHashSplit->{ $templateId->{$k}->[0] }->{ $msr[0] }  = [$k, $templateId->{$k}->[0]];

                    #$eachTemplateHash->{ $templateId->{$k}->[0] }->{ $eachTemplate->{modelreactions}->[$i]->{id}}  = [$k, $templateId->{$k}->[0]];
                    $eachTemplateHash->{ $templateId->{$k}->[0] }->{$msr[0]  }  = [$k, $templateId->{$k}->[0]];

                    if ($msr[0] eq 'rxn00000'){

                        $eachModelRxns->{$templateId->{$k}->[0]}++;

                        if ($eachTemplate->{modelreactions}->[$i]->{imported_gpr}){

                            $eachModelGARxns->{ $templateId->{$k}->[0] }->{$msr[0] } = [$eachTemplate->{modelreactions}->[$i]->{imported_gpr},$eachTemplate->{modelreactions}->[$i]->{name}];
                            $eachModelGARxnsCount->{$templateId->{$k}->[0]}++;


                        }
                        else{
                            push (@{$eachModelNonGARxns->{$templateId->{$k}->[0]}}, $msr[0] )

                        }


                    }
                    else {
                        #print "$eachTemplate->{modelreactions}->[$i]->{id}\t**$msr[0]**\t**$msr1[-1]**\n";
                        $eachModelRxns->{$templateId->{$k}->[0]}++;
                        $eachModelMSRxns->{$templateId->{$k}->[0]}++;

                        if ($eachTemplate->{modelreactions}->[$i]->{imported_gpr}){

                            $eachModelGARxns->{ $templateId->{$k}->[0] }->{$msr[0] } = [$eachTemplate->{modelreactions}->[$i]->{imported_gpr},$eachTemplate->{modelreactions}->[$i]->{name}];
                            $eachModelGARxnsCount->{$templateId->{$k}->[0]}++;
                            $eachModelMSGARxns->{ $templateId->{$k}->[0] }->{$msr[0]} = [$eachTemplate->{modelreactions}->[$i]->{imported_gpr},$eachTemplate->{modelreactions}->[$i]->{name}];
                            $eachModelMSGARCount->{$templateId->{$k}->[0]}++;

                        }
                        else{
                            push (@{$eachModelNonGARxns->{$templateId->{$k}->[0]}}, $msr[0] );

                        }

                    }

                } # for

                for (my $i=0; $i< @{$eachTemplate->{modelcompounds}}; $i++){
                    my @msc1 = split /\//, $eachTemplate->{modelcompounds}->[$i]->{compound_ref};
                    my @msc = split /_/, $msc1[-1];

                    if ($msc[0] eq 'cpd00000'){

                        $eachModelCpds->{$templateId->{$k}->[0]}++;
                    }
                    else {

                        $eachModelCpds->{$templateId->{$k}->[0]}++;
                        $eachModelMSCpds->{$templateId->{$k}->[0]}++;
                    }

                }

            };
            if ($@) {
              die "Error loading object from the workspace:\n".$@;
            }
        }

    }

    #print &Dumper ($eachModelRxns, $eachModelMSRxns, $eachModelGARxnsCount, $eachModelMSGARCount);


    my $userModelRxns;
    my $newModel;
    my $uMcounterHashGPR;
    my $uMcounterHashNOGPR;
    my $uMcounterHashGPRMS;
    my $uMcounterHashNOGPRNoMS;

    my $uMGPRrxns;
    my $uMNOGPRrxns;
    my $uMMSrxns;
    my $mMNonMSrxns;
    my $uMcounterHashGPRNoMS;
    my $uMNOGPRrxnsMS;

     my $userws = 'janakakbase:narrative_1518296212814';

 foreach my $k (keys $templateId){
    my $params->{output_model} = $k."_DM";

    print $params->{output_model} ."\n";

    eval {                                  #workspace=>$params->{workspace}, name=>$params->{output_model}
        $newModel = $wshandle->get_objects([{ workspace=>$userws, name=> $params->{output_model} }])->[0]{data};# ->{modelreactions};

        for (my $i=0; $i< @{$newModel->{modelreactions}}; $i++){
            $userModelRxns->{$newModel->{modelreactions}->[$i]->{id}} =1;
            push (@newModelArr,$newModel->{modelreactions}->[$i]->{id} );

            my @msr1 = split /\//, $newModel->{modelreactions}->[$i]->{reaction_ref};
            my @msr = split /_/, $msr1[-1];
            if ($msr[0] eq 'rxn00000'){
              $eachModelRxns->{ $params->{output_model}}++;
              #$mMNonMSrxns->{$params->{output_model}}->{$newModel->{modelreactions}->[$i]->{id}} =1;
              $mMNonMSrxns->{$params->{output_model}}->{$msr[0]} =1;
            }
            else {

              $eachModelRxns->{ $params->{output_model}}++;
              $eachModelMSRxns->{ $params->{output_model}}++;
              $uMMSrxns->{$params->{output_model}}->{$msr[0]} =1;
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

}

    #print &Dumper ($eachModelRxns, $eachModelMSRxns, $eachModelGARxnsCount, $eachModelMSGARCount);




    my $uMcounterHashNOGPRMS;
    my $uMNOGPRrxnsMS;
    my $uMcounterHashNOGPRNoMS;
    my $uMNOGPRrxnsNoMS;
    my $uMcounterHashGPRMS;
    my $uMGPRrxnsMS;
    my $uMGPRrxnsNoMS;

        #####
        foreach my $k (keys $templateId){
            my $params->{output_model} = $k."_DM";
            my $um = $k."_DM";
            print $params->{output_model} ."\n";
            $newModel = $wshandle->get_objects([{ workspace=>$userws, name=> $um }])->[0]{data};# ->{modelreactions};

            my $uMgprRxnCount =0;
            for (my $i=0; $i< @{$newModel->{modelreactions}}; $i++){
                my @msr1 = split /\//, $newModel->{modelreactions}->[$i]->{reaction_ref};
                my @msr = split /_/, $msr1[-1];

                my $gprCheck = $newModel->{modelreactions}->[$i];
                my @rid = split /_/, $gprCheck->{id};
                if (!@{$gprCheck->{modelReactionProteins} }){

                        if (exists $uMMSrxns->{$um}->{ $msr[0] }){
                            $uMcounterHashNOGPRMS->{$um}++;
                            $uMNOGPRrxnsMS->{$msr[0]} = 1;
                            #print "NO GPR $gprCheck->{id}\t $um \n";
                        }
                        elsif (exists $mMNonMSrxns->{ $um }->{ $msr[0]}){
                            $uMcounterHashNOGPRNoMS->{$um}++;
                            $uMNOGPRrxnsNoMS->{$msr[0]} = 1;
                            #print "NO GPR $gprCheck->{id}\t $um \n";
                        }
                        else{
                          next;
                        }



                }

                else{

                    my $gprFlag =0;
                    for (my $j=0; $j< @{$gprCheck->{modelReactionProteins}}; $j++){

                        my $subGPR = $gprCheck->{modelReactionProteins}->[$j];

                        for (my $n=0; $n< @{$subGPR->{modelReactionProteinSubunits}}; $n++){
                            my $eachGPR = $subGPR->{modelReactionProteinSubunits}->[$n];

                            if (!@{$eachGPR->{feature_refs} }){
                                next;

                            }
                            else{

                                $gprFlag =1;
                                last;
                            }
                        }
                    }
                    if ($gprFlag != 0){

                           if (exists $uMMSrxns->{$um}->{ $msr[0] }){
                                $uMcounterHashGPRMS->{$um}++;
                                $uMGPRrxnsMS->{$um}->{$msr[0]} = 1;
                                $uMgprRxnCount++;
                                #print "GPR $msr[0]\t $um\t$templateId->{$k}->[1] \n";
                            }
                            elsif (exists $mMNonMSrxns->{ $um }->{ $msr[0]}){
                                $uMcounterHashGPRNoMS->{$um}++;
                                $uMGPRrxnsNoMS->{$msr[0]} = 1;
                                $uMgprRxnCount++
                                #print "GPR $gprCheck->{id}\t $templateId->{$k}->[1] \n";
                            }
                            else{
                              next;
                            }

                    }
                    else{

                            if (exists $uMMSrxns->{$um}->{ $msr[0] }){
                                $uMcounterHashNOGPRMS->{$um}++;

                                $uMNOGPRrxnsMS->{$msr[0]} = 1;
                                #print "NO GPR $gprCheck->{id}\t $um \n";
                            }
                            elsif (exists $mMNonMSrxns->{ $um }->{ $msr[0]}){
                                $uMcounterHashNOGPRNoMS->{$um}++;

                                $uMNOGPRrxnsNoMS->{$msr[0]} = 1;
                                #print "NO GPR $gprCheck->{id}\t $um \n";
                            }
                            else{
                              next;
                            }

                    }
                }#else
            }#for

        }


 #print &Dumper ($uMcounterHashNOGPRMS, $uMcounterHashNOGPRNoMS);

 #print &Dumper ($eachModelRxns, $eachModelMSRxns, $eachModelGARxnsCount, $eachModelMSGARCount, $uMcounterHashGPRMS, $uMcounterHashGPRNoMS, $uMcounterHashNOGPRNoMS );
#header line
print "Model Name\tTotal Rxns\tTotal MS Rxns\tGARxnCount(templateModel)\tGARxnCount-MSEED (templateModel)\tGARxncount-MSEED\tGARxncount-NonMSEED\tRxns-NOGPR-NO_MSSEED\n";
foreach my $k (sort keys $templateId){
    my $um = $k."_DM";
    my $tm = $k."_KBase";

    #template model stats
    #print "$tm\t$eachModelRxns->{$tm}\t$eachModelMSRxns->{$tm}\t$eachModelGARxnsCount->{$tm}\t$eachModelMSGARCount->{$tm}\t$uMcounterHashGPRMS->{$tm}\t$uMcounterHashGPRNoMS->{$tm}\t$uMcounterHashNOGPRNoMS->{$tm}\n";

    #user model stats
    print "$um\t$eachModelRxns->{$um}\t$eachModelMSRxns->{$um}\t$eachModelGARxnsCount->{$um}\t$eachModelMSGARCount->{$um}\t$uMcounterHashGPRMS->{$um}\t$uMcounterHashGPRNoMS->{$um}\t$uMcounterHashNOGPRNoMS->{$um}\n";

}


my $fractionUMGPRMSRxnsFromTemplate;



foreach my $k (sort keys $eachModelMSGARxns){

    #$uMGPRrxnsMS
    my $um = $k;
    $um  =~ s/_KBase/_DM/g;
    foreach my $r (keys $eachModelMSGARxns->{$k}){

        #print "$um\t$k\t$r\n";

        if (exists $uMGPRrxnsMS->{$um}->{$r}){

            $fractionUMGPRMSRxnsFromTemplate->{$um}++;

            #print "$um\t$r\n";


        }

    }

}

#print &Dumper ($fractionUMGPRMSRxnsFromTemplate);



die;

    #END build_model_stats
    my @_bad_returns;
    (ref($output) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"output\" (value was \"$output\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to build_model_stats:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'build_model_stats');
    }
    return($output);
}




=head2 update_model

  $output = $obj->update_model($params)

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

sub update_model
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to update_model:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'update_model');
    }

    my $ctx = $kb_fungalmodeling::kb_fungalmodelingServer::CallContext;
    my($output);
    #BEGIN update_model
    my $token=$ctx->token;
    my $provenance=$ctx->provenance;
    my $output = "stats_return";
    print("Starting updating models:\n");
    print(Dumper($params) . "\n");


    my $fbaO = new fba_tools::fba_toolsClient( $self->{'callbackURL'},
                                                            ( 'service_version' => 'release',
                                                              'async_version' => 'release',
                                                            )
                                                           );



    my $wshandle= Workspace::WorkspaceClient->new($self->{'workspace-url'},token=>$token);


    my $template_ws = 'jplfaria:narrative_1510597445008';# 'janakakbase:narrative_1513399583946'; # template workspaces
    my $template_genome_ref = 'FungalTemplate.genome';
    my $template_model_ref = 'master_fungal_template';
    my $ws_name = $params->{workspace};
    my $protCompId = 'proteinComp'.$params->{genome_ref};
    my $tmpGenome;
    my $genome_list;

    #production ws = 28129

    #ci ws = 32690 Fungal_Genomes
    #ci ws test workspace with the Narrative = 29112

      eval{
        #$genome_list = $wshandle->list_objects({ids=>[28129],type=>"KBaseGenomes.Genome"});
        $genome_list = $wshandle->list_objects({ids=>[32690],type=>"KBaseGenomes.Genome"});

      };
      if ($@) {
        die "Error loading fungal genome list from workspace:\n".$@;
      }


      open INFILE, "/kb/module/data/genomes_species_names_taxonomy_info.txt" or die "Couldn't open meta file $!\n";

      my %metaHash;
      while (defined(my $inputMeta = <INFILE>)){
        chomp $inputMeta;
        my @orgName = split /\t/, $inputMeta;
        #print "*$input*\n";
        $metaHash{$orgName[0]} = [$orgName[1], $orgName[2]];
      }
      close INFILE;

      open INFILETax, "/kb/module/data/taxa_information_modified" or die "Couldn't open taxonomy file $!\n";

      my %taxaHash;
      while (defined(my $inputMeta = <INFILETax>)){
        chomp $inputMeta;
        my @orgName = split /,/, $inputMeta;
        #print "*$input*\n";
        #if (defined $orgName[4]){
            $taxaHash{$orgName[0]} = [$orgName[1], $orgName[3], $orgName[4], $orgName[2] ];
        #}

      }
      close INFILETax;
      print &Dumper (\%taxaHash);



      #This section refers to change the taxonomy
      foreach my $p (@$genome_list){


          #my $eachGenome = $wshandle->get_objects([{workspace=>$p->[7],name=>$p->[1]}] )->[0];
          #my $pubGenomeRef = $p->[6]."/".$p->[0]."/".$p->[4];

          #$eachGenome->{info}->[10]->{Name} = 'Xylona heveae TC161';
          #print &Dumper ($eachGenome->{info});
          #print $p->[7];

          if (exists $taxaHash{$p->[1]} ){ #&& exists $metaHash{$p->[1]}){

            my $eachGenome = $wshandle->get_objects([{workspace=>$p->[7],name=>$p->[1]}] )->[0];
          #my $pubGenomeRef = $p->[6]."/".$p->[0]."/".$p->[4];

          #$eachGenome->{info}->[10]->{Name} = 'Xylona heveae TC161';
            print &Dumper ($eachGenome->{info});

            print "\n$p->[1]\t $taxaHash{$p->[1]}->[0]\t$taxaHash{$p->[1]}->[1]\t$taxaHash{$p->[1]}->[2]\n";
            $eachGenome->{data}->{taxonomy} =  $taxaHash{$p->[1]}->[1].", ".$taxaHash{$p->[1]}->[2];
            $eachGenome->{data}->{scientific_name} =  $taxaHash{$p->[1]}->[0];
            $eachGenome->{data}->{source} = "MycoCosm";
            $eachGenome->{data}->{domain} = "Eukaryota";
            $eachGenome->{data}->{source_id} = $p->[1];



              my $obj_info_list = undef;
              eval {
                  $obj_info_list = $wshandle->save_objects({
                    'workspace'=>$p->[7],
                    'objects'=>[{
                    'type'=>'KBaseGenomes.Genome',
                    'data'=>$eachGenome->{data},
                    'name'=> $p->[1],
                    'provenance'=>$provenance
                  }]
                  });
              };
              if ($@) {
                  die "Error saving modified genome object to workspace:\n".$@;
              }

              print &Dumper ($obj_info_list);

          }

          else {

              print "Genome  $p->[1] was not present\n";

          }

      }

=head # This sectioin refers to change the genome name
      foreach my $p (@$genome_list){

          my $eachGenome = $wshandle->get_objects([{workspace=>$p->[7],name=>$p->[1]}] )->[0];
          #my $pubGenomeRef = $p->[6]."/".$p->[0]."/".$p->[4];

          #$eachGenome->{info}->[10]->{Name} = 'Xylona heveae TC161';
          print &Dumper ($eachGenome->{info});
          #print $p->[7];

          if (exists $metaHash{$p->[1]} ){

            print "\n$p->[1]\t $metaHash{$p->[1]}->[0]\n";
            $eachGenome->{data}->{scientific_name} =  $metaHash{$p->[1]}->[0];



              my $obj_info_list = undef;
              eval {
                  $obj_info_list = $wshandle->save_objects({
                    'workspace'=>$p->[7],
                    'objects'=>[{
                    'type'=>'KBaseGenomes.Genome',
                    #'meta'=>$eachGenome->{info}->[10],
                    'data'=>$eachGenome->{data},
                    'name'=> $p->[1],
                    'provenance'=>$provenance
                  }]
                  });
              };
              if ($@) {
                  die "Error saving modified genome object to workspace:\n".$@;
              }

              print &Dumper ($obj_info_list);

          }

          else {

              print "Genome  $p->[1] was not present\n";

          }

      }
=cut

    #END update_model
    my @_bad_returns;
    (ref($output) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"output\" (value was \"$output\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to update_model:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'update_model');
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
media_ref has a value which is a string
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
media_ref has a value which is a string
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
