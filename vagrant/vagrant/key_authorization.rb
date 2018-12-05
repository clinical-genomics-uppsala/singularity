def authorize_key(config, user="root", *key_paths)
  [*key_paths, nil].each do |key_path|
    if key_path.nil?
      fail "Public key not found at following paths: #{key_paths.join(', ')}"
    end

    full_key_path = File.expand_path(key_path)

    if File.exists?(full_key_path)
      config.vm.provision 'file',
        run: 'once',
        source: full_key_path,
        destination: '/tmp/root_pubkey'

     puts "User is #{user}"

     if user == "root" then
      config.vm.provision 'shell',
        privileged: true,
        run: 'once',
        inline:
          "echo \"Creating /root/.ssh/authorized_keys with #{key_path}\" && " +
	  "mkdir -p /root/.ssh/ && " +
          'rm -f /root/.ssh/authorized_keys && ' +
          'mv /tmp/root_pubkey /root/.ssh/authorized_keys && ' +
          'chown root:root /root/.ssh/authorized_keys && ' +
          'chmod 600 /root/.ssh/authorized_keys && ' +
          'rm -f /tmp/root_pubkey && ' +
          'echo "Done!"'
     else
      config.vm.provision 'shell',
        privileged: true,
        run: 'once',
        inline:
          "echo \"Creating /home/#{user}/.ssh/authorized_keys with #{key_path}\" && " +
	  "mkdir -p /home/#{user}/.ssh/ && " +
          "rm -f /home/#{user}/.ssh/authorized_keys && " +
          "mv /tmp/root_pubkey /home/#{user}/.ssh/authorized_keys && " +
          "chown #{user}:#{user} /home/#{user}/.ssh/authorized_keys && " +
          "chmod 600 /home/#{user}/.ssh/authorized_keys && " +
          "rm -f /tmp/root_pubkey && " +
          'echo "Done!"'
     end

     break
    end
  end
end
