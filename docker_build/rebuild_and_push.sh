docker login
docker build -t eagle .
docker tag eagle loven7doo/eagle
docker push loven7doo/eagle
