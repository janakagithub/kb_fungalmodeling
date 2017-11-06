use kb_fungalmodeling::kb_fungalmodelingImpl;

use kb_fungalmodeling::kb_fungalmodelingServer;
use Plack::Middleware::CrossOrigin;



my @dispatch;

{
    my $obj = kb_fungalmodeling::kb_fungalmodelingImpl->new;
    push(@dispatch, 'kb_fungalmodeling' => $obj);
}


my $server = kb_fungalmodeling::kb_fungalmodelingServer->new(instance_dispatch => { @dispatch },
				allow_get => 0,
			       );

my $handler = sub { $server->handle_input(@_) };

$handler = Plack::Middleware::CrossOrigin->wrap( $handler, origins => "*", headers => "*");
