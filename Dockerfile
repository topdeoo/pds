from ubuntu:lastest
run sed -i 's@//.*archive.ubuntu.com@//mirrors.ustc.edu.cn@g' /etc/apt/sources.list \
    && sed -i 's/security.ubuntu.com/mirrors.ustc.edu.cn/g' /etc/apt/sources.list \
    && sed -i 's/http:/https:/g' /etc/apt/sources.list \
    && apt-get update \
    && apt-get upgrade \
    && apt-get install -y neo-vim librange-v3-dev libboost-all-dev libtinyxml2-dev \
    && apt-get install -y curl wget tar git ca-certificates gpg gcc g++ make \
    && apt-get install -y clang clang-format clang-tidy clangd \

run wget https://packages.gurobi.com/10.0/gurobi10.0.1_linux64.tar.gz \
    && mkdir -p ~/opt/ \
    && mv gurobi10.0.1_linux64.tar.gz ~/opt/
    && cd ~/opt/ && tar -xzvf gurobi10.0.1_linux64.tar.gz \
    && rm gurobi10.0.1_linux64.tar.gz \
    && echo 'export GUROBI_HOME="/root/opt/gurobi100/linux64"\n \
            export GRB_LICENSE_FILE="/root/opt/gurobi.lic"\n \
            export PATH="${PATH}:${GUROBI_HOME}/bin" \n \
            export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"\n' >> ~/.bashrc \
    && source ~/.bashrc

run wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null \
    && echo 'deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ jammy main' | tee /etc/apt/sources.list.d/kitware.list >/dev/null \
    && apt-get update && rm /usr/share/keyrings/kitware-archive-keyring.gpg \
    && apt-get install kitware-archive-keyring \
    && apt-get install cmake

copy ~/pds /root/pds
