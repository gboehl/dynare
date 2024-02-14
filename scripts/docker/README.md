# Dynare Docker Containers
We provide a range of pre-configured Docker containers for [Dynare](https://dynare.org), which include both Octave and MATLAB (pre-configured with Dynare already in the PATH) and all recommended toolboxes. These containers are ideal for using Dynare in CI/CD environments ([example Workflow for GitHub Actions](https://github.com/JohannesPfeifer/DSGE_mod/tree/master/.github/workflows)) or High Performance Computing clusters with either [Docker, ENROOT or Singularity](https://wiki.bwhpc.de/e/BwUniCluster2.0/Containers).

To minimize maintenance efforts while ensuring high levels of security, reliability, and performance, our Docker containers are built from the official [MATLAB container base image](https://hub.docker.com/r/mathworks/matlab) using a custom [Dockerfile](Dockerfile). For more information on building and customizing the containers, see the [built instructions and customization](#built-instructions-and-customization) section below. Additionally, we provide an example [docker-compose file](docker-compose.yml) for complete access to the Ubuntu Desktop via VNC.

## Supported tags

| Tags   | Dynare Version | MATLAB® Version | Octave Version | Operating System | Base Image              |
|--------|----------------|-----------------|----------------|------------------|-------------------------|
| latest | 6.0            | R2023b          | 8.4.0 (PPA)    | Ubuntu 22.04     | mathworks/matlab:R2023b |
| 6.0    | 6.0            | R2023b          | 8.4.0 (PPA)    | Ubuntu 22.04     | mathworks/matlab:R2023b |
| 5.5    | 5.5            | R2023b          | 6.4.0          | Ubuntu 22.04     | mathworks/matlab:R2023b |
| 5.4    | 5.4            | R2023a          | 5.2.0          | Ubuntu 20.04     | mathworks/matlab:R2023a |
| 5.3    | 5.3            | R2022b          | 5.2.0          | Ubuntu 20.04     | mathworks/matlab:R2022b |
| 5.2    | 5.2            | R2022a          | 5.2.0          | Ubuntu 20.04     | mathworks/matlab:R2022a |
| 5.1    | 5.1            | R2022a          | 5.2.0          | Ubuntu 20.04     | mathworks/matlab:R2022a |
| 5.0    | 5.0            | R2021b          | 5.2.0          | Ubuntu 20.04     | mathworks/matlab:R2021b |
| 4.6.4  | 4.6.4          | R2021a          | 5.2.0          | Ubuntu 20.04     | mathworks/matlab:R2021a |

Note that we use an inofficial [PPA](https://launchpad.net/~ubuntuhandbook1/+archive/ubuntu/octave) (maintained by [https://ubuntuhandbook.org](https://ubuntuhandbook.org)) to install Octave 8.4.0 on Ubuntu 22.04, the usual disclaimer on PPAs applies.
Once Ubuntu 24.04 is released, we will switch to the version from the official repositories.

## How to interact with the container

To pull the latest image to your machine, execute:
```sh
docker pull dynare/dynare:latest
```
or a specific version:
```sh
docker pull dynare/dynare:6.0
```
In the following we assume that you have access to a MATLAB license and show different workflows how to interact with the container.
Obviously, you need to adjust the environmental variable `MLM_LICENSE_FILE` to your use-case, please refer to the [MATLAB license](#matlab-license) section on licensing information and how to pass a personal license.
Alternatively, if you don't have access to a license or the closed-source mentality of MATLAB is an irreconcilable issue for you, you can equally well use Dynare with the free and open-source alternative Octave.

### Run Dynare in an interactive MATLAB session in the browser

To launch the container with the `-browser` option, execute:
```sh
docker run -it --rm -p 8888:8888 -e MLM_LICENSE_FILE=27000@matlab-campus.uni-tuebingen.de --shm-size=512M dynare/dynare:latest -browser
```
You will receive a URL to access MATLAB in a web browser, for example: `http://localhost:8888` or another IP address that you can use to reach your server, such as through a VPN like [Tailscale](https://tailscale.com) if you are behind a firewall. Enter the URL provided into a web browser. Note that if you set `MLM_LICENSE_FILE` to empty or leave it out from the command, you will be prompted to enter credentials for a MathWorks account associated with a MATLAB license. If you are using a network license manager, switch to the Network License Manager tab and enter the license server address instead. After providing your license information, a MATLAB session will start in the browser. This may take several minutes. To modify the behavior of MATLAB when launched with the `-browser` flag, pass environment variables to the `docker run` command. For more information, see [Advanced Usage](https://github.com/mathworks/matlab-proxy/blob/main/Advanced-Usage.md).

Note that the `-browser` flag is supported by base images starting from `mathworks/matlab:R2022a` using [noVNC](https://novnc.com). Some browsers may not support this workflow.

### Run Ubuntu desktop and interact with it via VNC

To start the Ubuntu desktop with activated VNC server, execute:
```sh
docker run -it --rm -p 5901:5901 -p 6080:6080 -e PASSWORD=dynare -e MLM_LICENSE_FILE=27000@matlab-campus.uni-tuebingen.de --shm-size=512M dynare/dynare:latest -vnc
```
To connect to the Ubuntu desktop, either:

- Point a browser to port 6080 of the docker host machine running this container (`http://hostname:6080`).
- Use a VNC client to connect to display 1 of the docker host machine (`hostname:1`). The VNC password is `matlab` by default, you can change that by adjusting the `PASSWORD` environment variable in the run command.
- If you are behind a firewall, we recommend to use a VPN such as [Tailscale](https://tailscale.com) such that you can access the VNC server via the Tailscale address of the server. 

### Run Dynare with Octave in an interactive command prompt

To start the container and run Dynare with Octave in an interactive command prompt, execute:
```sh
docker run -it --rm --shm-size=512M dynare/dynare:latest octave
```

### Run Dynare with MATLAB non-interactively in batch mode

To start the container and run an example mod file using Dynare with MATLAB execute:
```sh
docker run --rm -e MLM_LICENSE_FILE=27000@matlab-campus.uni-tuebingen.de dynare/dynare:latest matlab -batch "cd dynare/examples; dynare example1"
```

### Run a bash shell inside the container

To start a bash shell inside the container, execute:

```sh
docker run -it --rm --shm-size=512M dynare/dynare:latest -shell
```
You can also non-interactively run a sequence of commands:
```sh
docker run --rm --shm-size=512M \
  -e MLM_LICENSE_FILE=27000@matlab-campus.uni-tuebingen.de \
  dynare/dynare:latest /bin/bash -c "\
    cd /home/matlab/dynare/examples && \
    matlab -batch 'dynare example1 console' && \
    octave --eval 'dynare example1 console'"
```

### Run MATLAB Desktop using X11

To start the container and run MATLAB Desktop using X11, execute:
```sh
xhost +
docker run -it --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix:ro -e MLM_LICENSE_FILE=27000@matlab-campus.uni-tuebingen.de --shm-size=512M dynare/dynare:latest
```
The Desktop window of MATLAB will open on your machine. Note that the command above works only on a Linux operating system with X11 and its dependencies installed.

## Additional information

### MATLAB license

To run this container, your license must be [configured for cloud use](https://mathworks.com/help/install/license/licensing-for-mathworks-products-running-on-the-cloud.html). Individual and Campus-Wide licenses are already configured for cloud use. If you have a different license type, please contact your license administrator to configure it for cloud use. You can identify your license type and administrator by viewing your MathWorks Account. Administrators can consult the "Administer Network Licenses" documentation. If you don't have a MATLAB license, you can obtain a trial license at [MATLAB Trial for Docker](https://de.mathworks.com/campaigns/products/trials/targeted/dkr.html).
Lastly, if you run the container via a GitHub workflow, you don't need to provide a license as the IP of the GitHub runner is already covered by a sponsored MATLAB license.

#### Network license
If you're using a network license, you can pass the port and hostname via the `MLM_LICENSE_FILE` environmental variable in your `docker run` command or Docker Compose file. Here's an example `docker run` command that uses a network license:
```sh
docker run --init -it --rm -e MLM_LICENSE_FILE=27000@matlab-campus.uni-tuebingen.de --shm-size=512M dynare/dynare:latest matlab -batch "cd dynare/examples; dynare example1"
```

#### Personal License

To use a personal license, you must first create a license file via the MATHWORKS License Center, refer to [Option 2](https://de.mathworks.com/matlabcentral/answers/235126-how-do-i-generate-a-matlab-license-file#answer_190013) for detailed instructions.
For this process, you will need the `username` and a `host ID`. In the container, the username is predefined as `matlab`.
The `host ID` corresponds to the MAC address of any network adapter in the container.
In Docker, you can supply a [randomly generated MAC address](https://miniwebtool.com/mac-address-generator/) (e.g., A6-7E-1A-F4-9A-92) during the docker run command.
Download the file from MATHWORKS License Center and ensure you provide the container with access to the license file by mounting it as a (read-only) volume.

Here is an example `docker run` command that utilizes a license file named `license.lic`, which is located in a folder `$HOME/matlab-license` on the host machine; the MAC address associated with the license is set to `A6-7E-1A-F4-9A-92`:
```sh
docker run --init -it --rm --mac-address A6-7E-1A-F4-9A-92 --shm-size=512M -v $HOME/matlab-license/:/licenses:ro -e MLM_LICENSE_FILE=/licenses/license.lic dynare/dynare:latest matlab -batch "cd dynare/examples; dynare example1"
```

### Environment variables

When running the `docker run` command, you can specify environment variables using the `-e` option. The [base image](https://hub.docker.com/r/mathworks/matlab) documentation lists the available variables, such as `MLM_LICENSE_FILE`, `PASSWORD`, or `PROXY_SETTINGS`.

## Built instructions and customization

Here are the commands to create the Docker images available at [Docker Hub](https://hub.docker.com/r/dynare/dynare):
```sh
docker build --build-arg MATLAB_RELEASE=R2023b --build-arg DYNARE_RELEASE=6.0 -t dynare/dynare:latest .
docker build --build-arg MATLAB_RELEASE=R2023b --build-arg DYNARE_RELEASE=6.0 -t dynare/dynare:6.0 .
docker build --build-arg MATLAB_RELEASE=R2023b --build-arg DYNARE_RELEASE=5.5 -t dynare/dynare:5.5 .
docker build --build-arg MATLAB_RELEASE=R2023a --build-arg DYNARE_RELEASE=5.4 -t dynare/dynare:5.4 .
docker build --build-arg MATLAB_RELEASE=R2022b --build-arg DYNARE_RELEASE=5.3 -t dynare/dynare:5.3 .
docker build --build-arg MATLAB_RELEASE=R2022a --build-arg DYNARE_RELEASE=5.2 -t dynare/dynare:5.2 .
docker build --build-arg MATLAB_RELEASE=R2022a --build-arg DYNARE_RELEASE=5.1 -t dynare/dynare:5.1 .
docker build --build-arg MATLAB_RELEASE=R2021b --build-arg DYNARE_RELEASE=5.0 -t dynare/dynare:5.0 .
docker build --build-arg MATLAB_RELEASE=R2021a --build-arg DYNARE_RELEASE=4.6.4 -t dynare/dynare:4.6.4 .
```

If you need to customize the container, there are two ways to do so. You can either adjust the [Dockerfile](Dockerfile) and rebuild the container, or you can run the container interactively, make the necessary adjustments, and then commit the changes for later use. To commit changes to a container, use the `docker commit` command. This will create a new image with the changes you made. You can then use this image to start new containers with your customizations.
For more information on committing changes to a container, see the [Docker documentation](https://docs.docker.com/engine/reference/commandline/commit/) and how to [save changes in the containers](https://de.mathworks.com/help/cloudcenter/ug/save-changes-in-containers.html).

Note that if you plan to distribute your custom container, you should be aware of the licensing terms of any software included in the container.
The provided containers provide no inclusion or information about a MATLAB license file.

## License

This container includes commercial software products from The MathWorks, Inc. ("MathWorks Programs") and related materials. The MathWorks Programs are licensed under the MathWorks Software License Agreement, which is available in the MATLAB installation within this container.

The related materials in this container are licensed under separate licenses, which can be found in their respective folders.

Dynare is licensed under the GPL-3+.
