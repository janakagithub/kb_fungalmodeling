use strict;
use Data::Dumper;
use Test::More;
use Config::Simple;
use Time::HiRes qw(time);
use Bio::KBase::AuthToken;
use Workspace::WorkspaceClient;
use AssemblyUtil::AssemblyUtilClient;
use kb_fungalmodeling::kb_fungalmodelingImpl;

local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('kb_fungalmodeling');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = new Workspace::WorkspaceClient($ws_url,token => $token);
my $scratch = $config->{scratch};
my $callback_url = $ENV{'SDK_CALLBACK_URL'};
my $auth_token = Bio::KBase::AuthToken->new(token => $token, ignore_authrc => 1, auth_svc=>$config->{'auth-service-url'});
my $ctx = LocalCallContext->new($token, $auth_token->user_id);
$kb_fungalmodeling::kb_fungalmodelingServer::CallContext = $ctx;
my $impl = new kb_fungalmodeling::kb_fungalmodelingImpl();
=head
sub get_ws_name {
    if (!defined($ws_name)) {
        my $suffix = int(time * 1000);
        $ws_name = 'test_kb_fungalmodeling_' . $suffix;
        $ws_client->create_workspace({workspace => $ws_name});
    }
    return $ws_name;
}
=cut

#my $ws = "janakakbase:narrative_1509376805185";
#my $ws = "janakakbase:narrative_1498154949048";
my $ws = 'janakakbase:narrative_1509987427391';
my $input_genome = 'Aspergillus_terreus';
my $protInput = {
        workspace => $ws,
        genome_ref => $input_genome,
        template_model =>  'default_temp',
        translation_policy =>  'translate_only',
        gapfill_model => 0,
        output_model =>  'prpogated_model_out_'.$input_genome
};

my $template_ws= 'janakakbase:narrative_1509987427391';

my $templateBuild = {
        workspace => $template_ws,
        reference_genome => 'Neurospora_crassa',
        reference_model =>  'Neuropora_crassa_Model',
        output_model =>  'FungalTemplateModel'
};

eval {
    my $ret =$impl->build_fungal_model($protInput);
    #my $ret =$impl->build_fungal_template($templateBuild);
    print &Dumper ($ret);
};



my $err = undef;
if ($@) {
    $err = $@;
}
eval {
    if (defined($ws_name)) {
        $ws_client->delete_workspace({workspace => $ws_name});
        print("Test workspace was deleted\n");
    }
};
if (defined($err)) {
    use Scalar::Util 'blessed';
    if(blessed $err && $err->isa("Bio::KBase::Exceptions::KBaseException")) {
        die "Error while running tests. Remote error:\n" . $err->{data} .
            "Client-side error:\n" . $err;
    } else {
        die $err;
    }
}

{
    package LocalCallContext;
    use strict;
    sub new {
        my($class,$token,$user) = @_;
        my $self = {
            token => $token,
            user_id => $user
        };
        return bless $self, $class;
    }
    sub user_id {
        my($self) = @_;
        return $self->{user_id};
    }
    sub token {
        my($self) = @_;
        return $self->{token};
    }
    sub provenance {
        my($self) = @_;
        return [{'service' => 'kb_fungalmodeling', 'method' => 'please_never_use_it_in_production', 'method_params' => []}];
    }
    sub authenticated {
        return 1;
    }
    sub log_debug {
        my($self,$msg) = @_;
        print STDERR $msg."\n";
    }
    sub log_info {
        my($self,$msg) = @_;
        print STDERR $msg."\n";
    }
}
